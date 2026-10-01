#!/usr/bin/env python3
"""Retrieve published AlphaFold Database models from explicit protein identifiers.

This is a standalone, standard-library-only command. NCBI Gene IDs, RefSeq
proteins and INSDC (GenBank/ENA/DDBJ) protein accessions are resolved through
UniProt's ID mapping API. All mappings are retained without guessed ranking.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import os
import re
import shutil
import sys
import time
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urlparse
from urllib.request import HTTPRedirectHandler, Request, build_opener
from uuid import uuid4

LOGGER = logging.getLogger("afdb_downloader")
UNIPROT_HOST = "rest.uniprot.org"
AFDB_HOST = "alphafold.ebi.ac.uk"
ID_MAPPING_URL = f"https://{UNIPROT_HOST}/idmapping"
AFDB_API_URL = f"https://{AFDB_HOST}/api/prediction"
MAX_IDENTIFIERS = 500
MAX_MAPPINGS = 5000
MAX_REQUESTED_MODELS = 500
MAX_JSON_BYTES = 32 * 1024 * 1024
MAX_MODEL_BYTES = 100 * 1024 * 1024
MAX_PAE_BYTES = 200 * 1024 * 1024
MANIFEST_FIELDS = (
    "input_id",
    "input_type",
    "mapping_source",
    "uniprot_accession",
    "mapping_count",
    "status",
    "model_path",
    "named_model_path",
    "model_sha256",
    "model_bytes",
    "model_url",
    "entry_id",
    "model_version",
    "pae_status",
    "pae_path",
    "named_pae_path",
    "pae_sha256",
    "metadata_sequence_check",
    "detail",
)
_UNIPROT_PATTERN = re.compile(r"[A-Z][A-Z0-9]{5}(?:[A-Z0-9]{4})?")
_REFSEQ_PATTERN = re.compile(r"(?:NP|XP|YP|WP|ZP|AP)_\d+(?:\.\d+)?")
_INSDC_PROTEIN_PATTERN = re.compile(r"[A-Z]{3}(?:\d{5}|\d{7})(?:\.\d+)?")
_GENE_PATTERN = re.compile(r"[1-9]\d{0,11}")


class DownloadError(Exception):
    """Indicate invalid service data, unsupported identifiers, or failed I/O."""


class ModelUnavailable(DownloadError):
    """Indicate that AlphaFold Database has no canonical F1 model."""


@dataclass(frozen=True)
class ModelReference:
    """Describe a validated AlphaFold Database resource.

    Attributes:
        accession: UniProt accession assigned to the model.
        entry_id: Canonical AlphaFold Database F1 entry identifier.
        version: Model release version extracted from its filename.
        model_url: HTTPS URL for the selected coordinates.
        pae_url: HTTPS URL for PAE data, when available.
        sequence: UniProt sequence supplied by the model metadata, if present.
    """

    accession: str
    entry_id: str
    version: int
    model_url: str
    pae_url: str | None
    sequence: str | None


@dataclass(frozen=True)
class DownloadedAsset:
    """Describe a file written during this run.

    Attributes:
        path: Path relative to the output directory.
        digest: SHA-256 digest of the downloaded bytes.
        size: Number of bytes written.
    """

    path: str
    digest: str
    size: int


class _BoundedRedirects(HTTPRedirectHandler):
    """Reject redirects away from the two explicitly approved services."""

    def redirect_request(
        self,
        request: Request,
        file_pointer: Any,
        code: int,
        message: str,
        headers: Any,
        new_url: str,
    ) -> Request | None:
        """Validate the redirected URL before allowing urllib to follow it.

        Args:
            request: Original HTTP request.
            file_pointer: Response file object provided by urllib.
            code: HTTP redirect status code.
            message: HTTP status message.
            headers: Headers supplied with the redirect.
            new_url: Proposed redirect destination.

        Returns:
            A request for the approved redirect, or ``None`` if urllib
            declines the redirect.

        Raises:
            DownloadError: If the destination is outside the approved host.
        """
        original_host = urlparse(request.full_url).hostname
        validate_https_url(url=new_url, host=str(original_host))
        return super().redirect_request(
            request, file_pointer, code, message, headers, new_url
        )


def validate_https_url(*, url: str, host: str, suffix: str | None = None) -> str:
    """Accept an exact HTTPS host and, where needed, an exact file suffix.

    Args:
        url: URL supplied by a service.
        host: Only approved DNS host for this request.
        suffix: Optional expected file extension.

    Returns:
        The validated URL, unchanged.

    Raises:
        DownloadError: If the URL is not a permitted HTTPS resource.
    """
    parsed = urlparse(url)
    if (
        parsed.scheme != "https"
        or parsed.netloc != host
        or not parsed.path.startswith("/")
        or not parsed.path.strip("/")
        or parsed.params
        or parsed.query
        or parsed.fragment
        or (suffix is not None and not parsed.path.lower().endswith(suffix))
    ):
        raise DownloadError(f"Unapproved {host} resource URL: {url!r}")
    return url


class HttpClient:
    """Read bounded HTTPS resources with timeouts, retries, and safe redirects."""

    def __init__(self, *, timeout_seconds: float = 30.0, retries: int = 3) -> None:
        """Create a bounded client with an approved-host redirect handler.

        Args:
            timeout_seconds: Timeout for each HTTP attempt, in seconds.
            retries: Maximum retries following a transient HTTP failure.

        Raises:
            DownloadError: If timeout or retry settings are invalid.
        """
        if not math.isfinite(timeout_seconds) or timeout_seconds <= 0 or retries < 0:
            raise DownloadError("HTTP timeout must be positive; retries non-negative")
        self.timeout_seconds = timeout_seconds
        self.retries = retries
        self._opener = build_opener(_BoundedRedirects())

    def read(
        self,
        *,
        url: str,
        host: str,
        limit: int,
        form: dict[str, str] | None = None,
    ) -> bytes:
        """Read a bounded HTTPS response and retry transient failures.

        Args:
            url: HTTPS resource to request.
            host: Exact approved host for the request and any redirect.
            limit: Maximum number of response bytes to accept.
            form: Form fields for a POST request, or ``None`` for GET.

        Returns:
            The response bytes, without modification.

        Raises:
            ModelUnavailable: If the server reports HTTP 404.
            DownloadError: If the URL, response size, or request fails.
        """
        validate_https_url(url=url, host=host)
        body = None if form is None else urlencode(form).encode("ascii")
        for attempt in range(self.retries + 1):
            request = Request(
                url,
                data=body,
                headers={
                    "User-Agent": "afdb-ncbi-downloader/1.0 (research utility)",
                    "Accept": "application/json" if limit == MAX_JSON_BYTES else "*/*",
                },
                method="POST" if body is not None else "GET",
            )
            try:
                with self._opener.open(
                    request, timeout=self.timeout_seconds
                ) as response:
                    content = response.read(limit + 1)
                if len(content) > limit:
                    raise DownloadError(f"Response exceeds {limit} bytes: {url}")
                return content
            except HTTPError as error:
                if error.code == 404:
                    raise ModelUnavailable(f"HTTP 404: {url}") from error
                if error.code not in {429, 500, 502, 503, 504}:
                    raise DownloadError(f"HTTP {error.code} from {host}") from error
                reason = f"HTTP {error.code} from {host}"
            except (URLError, TimeoutError, OSError) as error:
                reason = f"Network error from {host}: {error}"
            if attempt == self.retries:
                raise DownloadError(f"{reason} after {attempt + 1} attempt(s)")
            delay = min(2**attempt, 8)
            LOGGER.warning("%s; retrying in %s second(s)", reason, delay)
            time.sleep(delay)
        raise AssertionError("Unreachable retry state")


def read_json(*, payload: bytes, label: str) -> dict[str, Any]:
    """Decode a UTF-8 JSON object returned by a remote service.

    Args:
        payload: Raw response bytes from the service.
        label: Service description used in validation errors.

    Returns:
        The decoded JSON object.

    Raises:
        DownloadError: If the bytes are not a valid JSON object.
    """
    try:
        document = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise DownloadError(f"Invalid JSON from {label}") from error
    if not isinstance(document, dict):
        raise DownloadError(f"Expected a JSON object from {label}")
    return document


def validate_identifier(*, value: str, identifier_type: str) -> str:
    """Validate and normalise an accession for its declared input type.

    Args:
        value: Submitted accession or numeric NCBI Gene ID.
        identifier_type: One of the supported ``--id-type`` choices.

    Returns:
        The stripped, uppercase identifier.

    Raises:
        DownloadError: If the identifier does not match its declared type.
    """
    normalised = value.strip().upper()
    patterns = {
        "uniprot": _UNIPROT_PATTERN,
        "refseq-protein": _REFSEQ_PATTERN,
        "genbank-protein": _INSDC_PROTEIN_PATTERN,
        "gene-id": _GENE_PATTERN,
    }
    if identifier_type == "ncbi-protein":
        valid = any(
            pattern.fullmatch(normalised) is not None
            for pattern in (_REFSEQ_PATTERN, _INSDC_PROTEIN_PATTERN)
        )
    else:
        valid = patterns[identifier_type].fullmatch(normalised) is not None
    if not valid:
        raise DownloadError(
            f"Invalid {identifier_type} identifier {value!r}; choose the "
            "correct --id-type and provide protein rather than nucleotide IDs"
        )
    return normalised


def read_identifiers(
    *, identifier_type: str, ids: Iterable[str], ids_file: Path | None
) -> tuple[str, ...]:
    """Read identifiers and keep unique entries in their input order.

    Args:
        identifier_type: Input type used to validate every identifier.
        ids: Identifiers supplied directly on the command line.
        ids_file: Optional UTF-8 file with one identifier per line.

    Returns:
        Unique, validated identifiers in their original order.

    Raises:
        DownloadError: If input cannot be read, is invalid, or exceeds
            the maximum number of identifiers.
    """
    requested = list(ids)
    if ids_file is not None:
        try:
            with ids_file.open(encoding="utf-8-sig") as stream:
                requested.extend(
                    line.strip()
                    for line in stream
                    if line.strip() and not line.lstrip().startswith("#")
                )
        except (OSError, UnicodeError) as error:
            raise DownloadError(
                f"Cannot read input file {ids_file}: {error}"
            ) from error
    if not requested:
        raise DownloadError("Supply --id or --ids-file")
    unique: list[str] = []
    seen: set[str] = set()
    for value in requested:
        identifier = validate_identifier(value=value, identifier_type=identifier_type)
        if identifier not in seen:
            unique.append(identifier)
            seen.add(identifier)
    if len(unique) > MAX_IDENTIFIERS:
        raise DownloadError(f"At most {MAX_IDENTIFIERS} distinct IDs per run")
    return tuple(unique)


def read_fasta(
    *, path: Path | None, identifiers: tuple[str, ...] | None
) -> dict[str, str]:
    """Read protein FASTA records keyed by their header accessions.

    Args:
        path: FASTA path, or ``None`` when no sequences were provided.
        identifiers: Required set of accessions, or ``None`` to accept all
            accessions from the FASTA file.

    Returns:
        Uppercase protein sequences keyed by uppercase accession; an empty
        mapping if ``path`` is ``None``.

    Raises:
        DownloadError: If the file is invalid, contains duplicate IDs, or
            does not match the requested accession set.
    """
    if path is None:
        return {}
    sequences: dict[str, str] = {}
    identifier: str | None = None
    parts: list[str] = []

    def store_sequence() -> None:
        """Validate the current FASTA record and store its sequence.

        Raises:
            DownloadError: If its sequence is invalid or ID is duplicated.
        """
        if identifier is None:
            return
        sequence = "".join(parts).upper().removesuffix("*")
        if not sequence or re.fullmatch(r"[A-Z]+", sequence) is None:
            raise DownloadError(f"Invalid protein FASTA sequence for {identifier}")
        if identifier in sequences:
            raise DownloadError(f"Duplicate FASTA identifier {identifier}")
        sequences[identifier] = sequence

    try:
        with path.open(encoding="utf-8-sig") as stream:
            for raw_line in stream:
                line = raw_line.strip()
                if line.startswith(">"):
                    store_sequence()
                    header = line[1:].split(maxsplit=1)
                    if not header:
                        raise DownloadError("Empty protein FASTA header")
                    identifier = header[0].upper()
                    parts = []
                elif line:
                    if identifier is None:
                        raise DownloadError("FASTA sequence appears before a header")
                    parts.append(line)
        store_sequence()
    except (OSError, UnicodeError) as error:
        raise DownloadError(f"Cannot read FASTA {path}: {error}") from error
    if identifiers is not None:
        missing = set(identifiers) - sequences.keys()
        unexpected = sequences.keys() - set(identifiers)
        if missing or unexpected:
            raise DownloadError(
                f"FASTA IDs must equal requested IDs; missing={sorted(missing)}, "
                f"unexpected={sorted(unexpected)}"
            )
    if not sequences:
        raise DownloadError("Protein FASTA contains no sequences")
    return sequences


def map_identifiers(
    *,
    identifiers: tuple[str, ...],
    identifier_type: str,
    client: HttpClient,
    wait_seconds: float = 180.0,
    polling_seconds: float = 2.0,
) -> dict[str, tuple[str, ...]]:
    """Map input identifiers to all their UniProt accessions.

    Args:
        identifiers: Validated input accessions in their original order.
        identifier_type: Declared type of the input identifiers.
        client: Bounded HTTP client used for UniProt mapping requests.
        wait_seconds: Maximum time allowed for each mapping job.
        polling_seconds: Delay between mapping status checks.

    Returns:
        UniProt accessions grouped by each input identifier. Unmapped
        identifiers have empty tuples.

    Raises:
        DownloadError: If an ID mapping job fails or returns invalid data.
    """
    if identifier_type == "uniprot":
        return {identifier: (identifier,) for identifier in identifiers}
    grouped: dict[str, list[str]] = {}
    for identifier in identifiers:
        grouped.setdefault(
            mapping_source_for(identifier=identifier, identifier_type=identifier_type),
            [],
        ).append(identifier)
    mapping: dict[str, tuple[str, ...]] = {}
    for source, group in grouped.items():
        LOGGER.info("Mapping %s input(s) using UniProt source %s", len(group), source)
        mapping.update(
            map_source_group(
                identifiers=tuple(group),
                source=source,
                client=client,
                wait_seconds=wait_seconds,
                polling_seconds=polling_seconds,
            )
        )
    return mapping


def mapping_source_for(*, identifier: str, identifier_type: str) -> str:
    """Select the UniProt mapping source for a validated input accession.

    Args:
        identifier: Validated accession or numeric NCBI Gene ID.
        identifier_type: Declared input type, including mixed NCBI proteins.

    Returns:
        UniProt mapping source name, or ``DIRECT_UNIPROT`` for UniProt IDs.

    Raises:
        DownloadError: If the declared input type is unsupported.
    """
    if identifier_type == "uniprot":
        return "DIRECT_UNIPROT"
    if identifier_type == "gene-id":
        return "GeneID"
    if identifier_type == "refseq-protein" or (
        identifier_type == "ncbi-protein"
        and _REFSEQ_PATTERN.fullmatch(identifier) is not None
    ):
        return "RefSeq_Protein"
    if identifier_type in {"genbank-protein", "ncbi-protein"}:
        return "EMBL-GenBank-DDBJ_CDS"
    raise DownloadError(f"Unsupported input type: {identifier_type}")


def map_source_group(
    *,
    identifiers: tuple[str, ...],
    source: str,
    client: HttpClient,
    wait_seconds: float,
    polling_seconds: float,
) -> dict[str, tuple[str, ...]]:
    """Run one asynchronous UniProt mapping job for a single ID source.

    Args:
        identifiers: Accessions accepted by the same mapping source.
        source: UniProt source database name for the mapping request.
        client: Bounded HTTP client used to submit and poll the job.
        wait_seconds: Maximum time to wait for the mapping result.
        polling_seconds: Delay between status requests.

    Returns:
        All unique UniProt matches per input accession, sorted by accession.

    Raises:
        DownloadError: If the job fails, times out, or returns invalid data.
    """
    submit = read_json(
        payload=client.read(
            url=f"{ID_MAPPING_URL}/run",
            host=UNIPROT_HOST,
            limit=MAX_JSON_BYTES,
            form={"from": source, "to": "UniProtKB", "ids": ",".join(identifiers)},
        ),
        label="UniProt ID mapping submission",
    )
    job_id = submit.get("jobId")
    if not isinstance(job_id, str) or re.fullmatch(r"[A-Za-z0-9]+", job_id) is None:
        raise DownloadError("UniProt did not return a valid ID mapping job ID")
    if (
        not math.isfinite(wait_seconds)
        or not math.isfinite(polling_seconds)
        or wait_seconds <= 0
        or polling_seconds <= 0
    ):
        raise DownloadError("Mapping time limits must be positive")
    deadline = time.monotonic() + wait_seconds
    while True:
        status = read_json(
            payload=client.read(
                url=f"{ID_MAPPING_URL}/status/{job_id}",
                host=UNIPROT_HOST,
                limit=MAX_JSON_BYTES,
            ),
            label="UniProt ID mapping status",
        )
        state = status.get("jobStatus")
        if state == "FAILED":
            raise DownloadError("UniProt ID mapping job failed")
        if state == "FINISHED" or "results" in status or "failedIds" in status:
            break
        if state not in {"RUNNING", "NEW", "PENDING"}:
            raise DownloadError(f"Unknown UniProt mapping state: {state!r}")
        if time.monotonic() + polling_seconds >= deadline:
            raise DownloadError("UniProt ID mapping timed out")
        time.sleep(polling_seconds)
    document = read_json(
        payload=client.read(
            url=f"{ID_MAPPING_URL}/stream/{job_id}",
            host=UNIPROT_HOST,
            limit=MAX_JSON_BYTES,
        ),
        label="UniProt ID mapping results",
    )
    results = document.get("results")
    if not isinstance(results, list) or len(results) > MAX_MAPPINGS:
        raise DownloadError("UniProt mapping results absent or above safety limit")
    matches: dict[str, set[str]] = {identifier: set() for identifier in identifiers}
    for row in results:
        if not isinstance(row, dict) or not isinstance(row.get("from"), str):
            raise DownloadError("Malformed UniProt mapping row")
        original = row["from"].upper()
        if original not in matches:
            raise DownloadError(f"Unexpected mapped input identifier: {original}")
        target = row.get("to")
        value = target.get("primaryAccession") if isinstance(target, dict) else target
        if not isinstance(value, str):
            raise DownloadError(f"Missing UniProt accession for {original}")
        accession = validate_identifier(value=value, identifier_type="uniprot")
        matches[original].add(accession)
    return {key: tuple(sorted(value)) for key, value in matches.items()}


def select_model(
    *, payload: bytes, accession: str, model_format: str
) -> ModelReference:
    """Select the latest unique canonical F1 model from AlphaFold metadata.

    Args:
        payload: Raw AlphaFold Database prediction metadata in JSON format.
        accession: UniProt accession expected in the metadata records.
        model_format: Coordinate format, either ``cif`` or ``pdb``.

    Returns:
        Validated reference to the latest matching canonical F1 model.

    Raises:
        ModelUnavailable: If no matching canonical model exists.
        DownloadError: If the metadata are invalid or latest records conflict.
    """
    try:
        records = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise DownloadError(f"Invalid AlphaFold metadata for {accession}") from error
    if not isinstance(records, list):
        raise DownloadError(f"AlphaFold metadata for {accession} is not a list")
    entry_id = f"AF-{accession}-F1"
    candidates: list[ModelReference] = []
    for record in records:
        if not isinstance(record, dict):
            raise DownloadError(f"Invalid AlphaFold entry for {accession}")
        if str(record.get("uniprotAccession", "")).upper() != accession:
            continue
        if record.get("entryId") != entry_id:
            continue
        field = "cifUrl" if model_format == "cif" else "pdbUrl"
        model_url = record.get(field)
        if not isinstance(model_url, str):
            continue
        validate_https_url(url=model_url, host=AFDB_HOST, suffix=f".{model_format}")
        filename = Path(urlparse(model_url).path).name
        pattern = rf"AF-{re.escape(accession)}-F1-model_v(\d+)\.{model_format}"
        match = re.fullmatch(pattern, filename)
        if match is None:
            raise DownloadError(f"Unexpected model filename for {accession}")
        version = int(match.group(1))
        raw_pae = record.get("paeDocUrl")
        pae_url = None
        if raw_pae:
            if not isinstance(raw_pae, str):
                raise DownloadError("Invalid PAE URL field")
            pae_url = validate_https_url(url=raw_pae, host=AFDB_HOST, suffix=".json")
        sequence = record.get("uniprotSequence")
        if sequence is not None and (
            not isinstance(sequence, str) or re.fullmatch(r"[A-Z]+", sequence) is None
        ):
            raise DownloadError("Invalid UniProt sequence in AlphaFold metadata")
        candidates.append(
            ModelReference(
                accession=accession,
                entry_id=entry_id,
                version=version,
                model_url=model_url,
                pae_url=pae_url,
                sequence=sequence,
            )
        )
    if not candidates:
        raise ModelUnavailable(f"No canonical F1 {model_format} for {accession}")
    latest = max(reference.version for reference in candidates)
    newest = [reference for reference in candidates if reference.version == latest]
    if any(reference != newest[0] for reference in newest[1:]):
        raise DownloadError(f"Conflicting latest AlphaFold records for {accession}")
    return newest[0]


def validate_model_bytes(*, payload: bytes, model_format: str) -> None:
    """Reject responses that are clearly not model coordinate files.

    Args:
        payload: Downloaded model response bytes.
        model_format: Expected coordinate format, ``cif`` or ``pdb``.

    Raises:
        DownloadError: If content is empty, binary, or lacks expected
            coordinate markers.
    """
    if not payload or b"\x00" in payload:
        raise DownloadError("Empty or binary coordinate response")
    if model_format == "cif":
        if not payload.lstrip().startswith(b"data_") or b"_atom_site." not in payload:
            raise DownloadError("Downloaded data are not a model mmCIF")
    elif not any(line.startswith(b"ATOM  ") for line in payload.splitlines()):
        raise DownloadError("Downloaded data are not a model PDB")


def write_asset(*, directory: Path, name: str, payload: bytes) -> DownloadedAsset:
    """Write a model asset atomically and calculate its SHA-256 digest.

    Args:
        directory: Existing canonical ``models`` output directory.
        name: Filename derived from a validated AlphaFold Database URL.
        payload: Downloaded model or PAE response bytes.

    Returns:
        Relative path, checksum, and byte size of the stored file.

    Raises:
        DownloadError: If the target file already exists.
        OSError: If the temporary file or atomic replacement cannot be written.
    """
    destination = directory / name
    if destination.exists():
        raise DownloadError(f"Output file already exists: {destination}")
    temporary = directory / f".{name}.{uuid4().hex}.tmp"
    try:
        with temporary.open("xb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return DownloadedAsset(
        path=f"models/{name}", digest=sha256(payload).hexdigest(), size=len(payload)
    )


def write_named_alias(
    *, output_dir: Path, asset: DownloadedAsset, identifier: str, accession: str
) -> str:
    """Expose a model under its submitted and UniProt IDs without re-fetching.

    Hard links avoid storing repeated large coordinate files; on file systems
    that do not permit links, the file is copied and published atomically.

    Args:
        output_dir: Root directory containing ``models`` and ``by_input``.
        asset: Canonical model or PAE asset already written to ``models``.
        identifier: Original validated input accession.
        accession: UniProt accession for the canonical asset.

    Returns:
        Path to the named file relative to ``output_dir``.

    Raises:
        DownloadError: If the named file already exists or cannot be written.
    """
    name = f"{identifier}__{accession}__{Path(asset.path).name}"
    destination = output_dir / "by_input" / name
    if destination.exists():
        raise DownloadError(f"Output file already exists: {destination}")
    temporary = destination.with_name(f".{name}.{uuid4().hex}.tmp")
    try:
        source = output_dir / asset.path
        try:
            os.link(source, temporary)
        except OSError:
            with source.open("rb") as original, temporary.open("xb") as copy:
                shutil.copyfileobj(original, copy)
                copy.flush()
                os.fsync(copy.fileno())
        os.replace(temporary, destination)
    except OSError as error:
        raise DownloadError(f"Could not name file for {identifier}: {error}") from error
    finally:
        temporary.unlink(missing_ok=True)
    return f"by_input/{name}"


def new_row(
    *, identifier: str, identifier_type: str, accession: str = "", count: int = 0
) -> dict[str, str]:
    """Create an initial manifest row for an input and optional mapping.

    Args:
        identifier: Validated original input identifier.
        identifier_type: Declared input identifier type.
        accession: Resolved UniProt accession, if available.
        count: Number of UniProt matches for the input identifier.

    Returns:
        Manifest row with every output field and default status values.

    Raises:
        DownloadError: If the identifier type has no mapping source.
    """
    row = {key: "" for key in MANIFEST_FIELDS}
    row.update(
        input_id=identifier,
        input_type=identifier_type,
        mapping_source=mapping_source_for(
            identifier=identifier, identifier_type=identifier_type
        ),
        uniprot_accession=accession,
        mapping_count=str(count),
        pae_status="NOT_REQUESTED",
        metadata_sequence_check="NOT_PROVIDED",
    )
    return row


def write_manifest(*, path: Path, rows: list[dict[str, str]]) -> None:
    """Write an auditable TSV manifest using an atomic file replacement.

    Tabs and newlines are removed from values to keep one output record per
    line, including when remote error messages contain control characters.

    Args:
        path: Output manifest path.
        rows: Model outcomes with keys from ``MANIFEST_FIELDS``.

    Raises:
        OSError: If the temporary manifest or replacement cannot be written.
    """
    temporary = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    try:
        with temporary.open("x", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(
                stream,
                fieldnames=MANIFEST_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for row in rows:
                writer.writerow(
                    {
                        key: re.sub(r"[\t\r\n]+", " ", str(row.get(key, "")))
                        for key in MANIFEST_FIELDS
                    }
                )
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def prepare_output(*, directory: Path) -> Path:
    """Create an empty output tree without overwriting previous results.

    Args:
        directory: Requested output directory, which must be new or empty.

    Returns:
        Resolved output directory containing ``models`` and ``by_input``.

    Raises:
        DownloadError: If the path is a symlink, a file, or a non-empty
            directory.
        OSError: If the output directories cannot be created.
    """
    expanded = directory.expanduser()
    if expanded.is_symlink():
        raise DownloadError("Output directory must not be a symlink")
    if expanded.exists():
        if not expanded.is_dir() or any(expanded.iterdir()):
            raise DownloadError("Output directory must be new or empty")
    else:
        expanded.mkdir(parents=True)
    (expanded / "models").mkdir()
    (expanded / "by_input").mkdir()
    return expanded.resolve()


def download_models(
    *,
    identifiers: tuple[str, ...],
    identifier_type: str,
    output_dir: Path,
    model_format: str,
    include_pae: bool,
    expected_sequences: dict[str, str],
    client: HttpClient,
    max_models: int = 50,
    mapper: Callable[..., dict[str, tuple[str, ...]]] = map_identifiers,
) -> list[dict[str, str]]:
    """Retrieve mapped models and write an outcome for every input mapping.

    Args:
        identifiers: Validated input identifiers in their original order.
        identifier_type: Declared type of those identifiers.
        output_dir: Prepared directory for models and the manifest.
        model_format: Requested coordinate format, ``cif`` or ``pdb``.
        include_pae: Whether to request available PAE JSON files.
        expected_sequences: Optional exact protein sequences by input ID.
        client: Bounded HTTP client for mapping and model requests.
        max_models: Maximum number of distinct UniProt models to retrieve.
        mapper: Callable that maps input IDs to UniProt accessions.

    Returns:
        Manifest rows, including unmapped, unavailable, and failed outcomes.
        The same rows are written to ``manifest.tsv``.

    Raises:
        OSError: If the manifest cannot be written.
    """
    try:
        if not 1 <= max_models <= MAX_REQUESTED_MODELS:
            raise DownloadError(
                f"--max-models must be between 1 and {MAX_REQUESTED_MODELS}"
            )
        mapping = mapper(
            identifiers=identifiers, identifier_type=identifier_type, client=client
        )
        distinct = {accession for matches in mapping.values() for accession in matches}
        if len(distinct) > max_models:
            raise DownloadError(
                f"Mapped {len(distinct)} distinct UniProt models, above the "
                f"--max-models limit of {max_models}; raise it explicitly or "
                "split the input list"
            )
    except DownloadError as error:
        rows = []
        for identifier in identifiers:
            row = new_row(identifier=identifier, identifier_type=identifier_type)
            row["status"] = "FAILED"
            row["detail"] = f"Identifier mapping failed: {error}"
            rows.append(row)
        write_manifest(path=output_dir / "manifest.tsv", rows=rows)
        return rows
    rows: list[dict[str, str]] = []
    references: dict[str, ModelReference | DownloadError] = {}
    downloaded: dict[str, tuple[DownloadedAsset, DownloadedAsset | None, str]] = {}
    for identifier in identifiers:
        accessions = mapping[identifier]
        if not accessions:
            row = new_row(identifier=identifier, identifier_type=identifier_type)
            row["status"] = "UNMAPPED"
            row["detail"] = "UniProt has no mapping for this identifier"
            rows.append(row)
            continue
        for accession in accessions:
            row = new_row(
                identifier=identifier,
                identifier_type=identifier_type,
                accession=accession,
                count=len(accessions),
            )
            if len(accessions) > 1:
                row["detail"] = "One of multiple UniProt mappings for this input"
            try:
                if accession not in references:
                    try:
                        payload = client.read(
                            url=f"{AFDB_API_URL}/{accession}",
                            host=AFDB_HOST,
                            limit=MAX_JSON_BYTES,
                        )
                        references[accession] = select_model(
                            payload=payload,
                            accession=accession,
                            model_format=model_format,
                        )
                    except DownloadError as error:
                        references[accession] = error
                reference = references[accession]
                if isinstance(reference, DownloadError):
                    raise reference
                row.update(
                    model_url=reference.model_url,
                    entry_id=reference.entry_id,
                    model_version=str(reference.version),
                )
                if identifier in expected_sequences:
                    if reference.sequence is None:
                        row["metadata_sequence_check"] = "UNAVAILABLE"
                        raise DownloadError(
                            "AlphaFold metadata has no sequence for comparison"
                        )
                    if expected_sequences[identifier] != reference.sequence:
                        row["metadata_sequence_check"] = "MISMATCH"
                        row["status"] = "SEQUENCE_MISMATCH"
                        row["detail"] += "; input FASTA differs from UniProt sequence"
                        rows.append(row)
                        continue
                    row["metadata_sequence_check"] = "MATCH"
                if accession not in downloaded:
                    model_bytes = client.read(
                        url=reference.model_url,
                        host=AFDB_HOST,
                        limit=MAX_MODEL_BYTES,
                    )
                    validate_model_bytes(payload=model_bytes, model_format=model_format)
                    filename = Path(urlparse(reference.model_url).path).name
                    model = write_asset(
                        directory=output_dir / "models",
                        name=filename,
                        payload=model_bytes,
                    )
                    pae = None
                    pae_status = "NOT_REQUESTED"
                    if include_pae:
                        pae_status = "NOT_AVAILABLE"
                        if reference.pae_url:
                            try:
                                pae_bytes = client.read(
                                    url=reference.pae_url,
                                    host=AFDB_HOST,
                                    limit=MAX_PAE_BYTES,
                                )
                                json.loads(pae_bytes.decode("utf-8"))
                                pae_name = Path(urlparse(reference.pae_url).path).name
                                pae = write_asset(
                                    directory=output_dir / "models",
                                    name=pae_name,
                                    payload=pae_bytes,
                                )
                                pae_status = "DOWNLOADED"
                            except (
                                DownloadError,
                                UnicodeError,
                                json.JSONDecodeError,
                            ) as error:
                                LOGGER.warning(
                                    "PAE unavailable for %s: %s", accession, error
                                )
                                pae_status = "FAILED"
                    downloaded[accession] = (model, pae, pae_status)
                model, pae, pae_status = downloaded[accession]
                named_model = write_named_alias(
                    output_dir=output_dir,
                    asset=model,
                    identifier=identifier,
                    accession=accession,
                )
                named_pae = (
                    write_named_alias(
                        output_dir=output_dir,
                        asset=pae,
                        identifier=identifier,
                        accession=accession,
                    )
                    if pae is not None
                    else ""
                )
                row.update(
                    status="DOWNLOADED",
                    model_path=model.path,
                    named_model_path=named_model,
                    model_sha256=model.digest,
                    model_bytes=str(model.size),
                    pae_status=pae_status,
                    pae_path="" if pae is None else pae.path,
                    named_pae_path=named_pae,
                    pae_sha256="" if pae is None else pae.digest,
                )
                LOGGER.info("Model %s retrieved for input %s", accession, identifier)
            except ModelUnavailable as error:
                row["status"] = "MODEL_NOT_AVAILABLE"
                row["detail"] += f"; {error}"
                LOGGER.info("Model missing for %s: %s", identifier, error)
            except DownloadError as error:
                row["status"] = "FAILED"
                row["detail"] += f"; {error}"
                LOGGER.error("Model retrieval failed for %s: %s", identifier, error)
            rows.append(row)
    write_manifest(path=output_dir / "manifest.tsv", rows=rows)
    return rows


def parse_args(*, argv: list[str] | None = None) -> argparse.Namespace:
    """Parse named command-line options for model retrieval.

    Args:
        argv: Explicit argument list, or ``None`` to read ``sys.argv``.

    Returns:
        Parsed command-line options.

    Raises:
        SystemExit: If argparse encounters invalid options or ``--help``.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--id-type",
        choices=(
            "uniprot",
            "refseq-protein",
            "genbank-protein",
            "ncbi-protein",
            "gene-id",
        ),
        required=True,
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--id", action="append", help="Repeat for several IDs")
    source.add_argument("--ids-file", type=Path, help="One ID per line, no header")
    source.add_argument(
        "--protein-fasta",
        type=Path,
        help="Protein FASTA with accession headers; verify sequences exactly",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--format", choices=("cif", "pdb"), default="cif")
    parser.add_argument("--include-pae", action="store_true")
    parser.add_argument(
        "--max-models",
        type=int,
        default=50,
        help=(
            "Maximum distinct models to retrieve "
            f"(default: 50; cap: {MAX_REQUESTED_MODELS})"
        ),
    )
    parser.add_argument(
        "--sequence-fasta",
        type=Path,
        help="Optional FASTA keyed by input IDs; exact metadata sequence required",
    )
    parser.add_argument("--timeout-seconds", type=float, default=30.0)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)


def main(*, argv: list[str] | None = None) -> int:
    """Run model retrieval and return a process exit status.

    Args:
        argv: Explicit argument list, or ``None`` to read ``sys.argv``.

    Returns:
        ``0`` for completed requests, ``1`` for recorded download or PAE
        failures, or ``2`` for invalid input or a fatal setup error.

    Raises:
        SystemExit: If argparse encounters invalid options or ``--help``.
    """
    args = parse_args(argv=argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    try:
        if args.protein_fasta is not None:
            if args.sequence_fasta is not None:
                raise DownloadError(
                    "Use --protein-fasta alone, or --ids-file with --sequence-fasta"
                )
            expected = read_fasta(path=args.protein_fasta, identifiers=None)
            identifiers = read_identifiers(
                identifier_type=args.id_type,
                ids=expected.keys(),
                ids_file=None,
            )
        else:
            identifiers = read_identifiers(
                identifier_type=args.id_type,
                ids=args.id or (),
                ids_file=args.ids_file,
            )
            expected = read_fasta(path=args.sequence_fasta, identifiers=identifiers)
        client = HttpClient(
            timeout_seconds=args.timeout_seconds,
            retries=args.retries,
        )
        output_dir = prepare_output(directory=args.output_dir)
        rows = download_models(
            identifiers=identifiers,
            identifier_type=args.id_type,
            output_dir=output_dir,
            model_format=args.format,
            include_pae=args.include_pae,
            expected_sequences=expected,
            client=client,
            max_models=args.max_models,
        )
    except (DownloadError, OSError) as error:
        LOGGER.error("Cannot complete request: %s", error)
        return 2
    counts = {
        status: sum(row["status"] == status for row in rows)
        for status in {row["status"] for row in rows}
    }
    LOGGER.info("Wrote %s; outcomes: %s", output_dir / "manifest.tsv", counts)
    return (
        1
        if counts.get("FAILED", 0) or any(row["pae_status"] == "FAILED" for row in rows)
        else 0
    )


if __name__ == "__main__":
    sys.exit(main())
