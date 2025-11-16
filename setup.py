from setuptools import setup
from setuptools.command.install import install as _install
import subprocess
import os
import sys

class InstallAndSetup(_install):
    def run(self):
        # Ubuntu-only guard
        os_release = '/etc/os-release'
        if not os.path.exists(os_release):
            raise SystemExit('Ubuntu only. /etc/os-release not found.')
        with open(os_release, 'r') as f:
            data = f.read()
        if 'ID=ubuntu' not in data:
            raise SystemExit('Ubuntu only. Detected non-Ubuntu OS.')

        _install.run(self)
        base_dir = os.path.abspath(os.path.dirname(__file__))
        sh = os.path.join(base_dir, '/src/setup.sh')
        if os.path.exists(sh):
            try:
                subprocess.check_call(['bash', sh], cwd=base_dir)
            except subprocess.CalledProcessError as e:
                raise SystemExit(f'setup.sh failed with exit code {e.returncode}')
        else:
            print('setup.sh not found, skipped running it.')

setup(
    name='alphagem-installer',
    version='1.0.0',
    description='Installer wrapper to prepare AlphaGEM tools and data (Ubuntu only).',
    long_description='Runs setup.sh after install to prepare external tools and data. This installer is intended for Ubuntu systems only.',
    long_description_content_type='text/plain',
    author='AlphaGEMs',
    packages=[],
    cmdclass={'install': InstallAndSetup},
    zip_safe=False,
)
