"""Maintenance entry point for the data cache: ``python -m lunarenv.data``.

The user-facing CLI (``lunarenv``) wraps the same :class:`DataStore` calls.
"""

from __future__ import annotations

import argparse
import sys

from lunarenv.data.store import DataStore


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m lunarenv.data", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    fetch = sub.add_parser("fetch", help="download catalog files or groups")
    fetch.add_argument("targets", nargs="+", help="catalog ids or group names")
    fetch.add_argument(
        "--pin", action="store_true", help="record checksums of unpinned files (maintainers)"
    )

    fetch.add_argument("--workers", type=int, default=4, help="files fetched at a time")
    fetch.add_argument(
        "--connections",
        type=int,
        default=1,
        help="range connections per file (workers x connections <= 8)",
    )

    sub.add_parser("verify", help="re-hash every cached file")
    sub.add_parser("list", help="list catalog groups and what is cached")

    args = parser.parse_args(argv)
    store = DataStore()

    if args.command == "fetch":
        for target in args.targets:
            if target in store.catalog:
                ids = [target]
            else:
                ids = [entry.id for entry in store.catalog.group(target)]
            records = store.fetch_many(
                ids, pin=args.pin, workers=args.workers, connections=args.connections
            )
            for record in records:
                print(f"{record.id}  {record.bytes:>12d}  {record.sha256}")
        return 0

    if args.command == "verify":
        status = store.verify()
        for entry_id, state in status.items():
            print(f"{entry_id}: {state}")
        return 0 if all(state in ("ok", "missing") for state in status.values()) else 1

    print(f"cache: {store.root}")
    for group in sorted(store.catalog.groups):
        entries = store.catalog.group(group)
        cached = sum(store.is_available(e.id) for e in entries)
        print(f"  {group:<26s} {cached}/{len(entries)} cached")
    return 0


if __name__ == "__main__":
    sys.exit(main())
