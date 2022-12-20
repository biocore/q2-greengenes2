from setuptools import setup, find_packages

import versioneer


install_requires = [
    "biom-format",
    "iow",
    "redbiom",
    "scikit-bio"
]


setup(
    name="q2-greengenes2",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    url="https://github.com/biocore/q2-greengenes2",
    license="BSD-3-Clause",
    description="Support methods for interaction with Greengenes2",
    entry_points={
        "qiime2.plugins":
        ["q2-gg2=q2_gg2.plugin_setup:plugin"]
    },
    package_data={
        'q2_gg2.tests': [],
        'q2_gg2': []
    },
    zip_safe=False,
    install_requires=install_requires,
)
