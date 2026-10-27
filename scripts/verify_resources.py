"""Verify the files supplied by the AlphaGEM runtime resource bundle."""

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def verify(root):
    root = Path(root).resolve()
    manifest = json.loads((root / 'resource-manifest.json').read_text())
    if manifest.get('schema_version') != 1:
        raise ValueError('Unsupported resource manifest version')
    seen = set()
    for row in manifest['files']:
        relative = Path(row['path'])
        path = (root / relative).resolve()
        if relative.is_absolute() or not path.is_relative_to(root) or row['path'] in seen:
            raise ValueError(f'Unsafe or duplicate resource path: {relative}')
        seen.add(row['path'])
        if not path.is_file() or path.stat().st_size != row['bytes'] or sha256(path) != row['sha256']:
            raise ValueError(f'Missing or changed resource: {relative}')
    if len(seen) != manifest['file_count']:
        raise ValueError('Resource file count does not match the manifest')
    return len(seen)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    try:
        count = verify(args.root)
    except (OSError, KeyError, ValueError) as error:
        parser.exit(1, f'Resource verification failed: {error}\n')
    print(f'Verified {count} resource files')


if __name__ == '__main__':
    main()
