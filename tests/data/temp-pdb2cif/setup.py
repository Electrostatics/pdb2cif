"""Set up the pdb2cif package."""
import setuptools


setuptools.setup(
    name="pdb2cif",
    version="0.0.1",
    description=(
        "This code converts to and from traditional PDB and new mmCIF formats"
    ),
    long_description=(
        "This code converts to and from the old PDB format "
        "(https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html) "
        "and the new mmCIF format (https://mmcif.wwpdb.org/)."
    ),
    python_requires=">=3.8",
    license="CC0-1.0",
    author="Nathan Baker",
    author_email="nathanandrewbaker@gmail.com",
    url="https://github.com/Electrostatics/pdb2cif",
    packages=setuptools.find_packages(),
    package_data={"": ["tests/data/*"]},
    install_requires=["mmcif_pdbx", "old_pdb>=1.1.0", "requests"],
    tests_require=["pytest"],
    keywords="science chemistry biophysics biochemistry",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Common Public License",
        "Natural Language :: English",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
