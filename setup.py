"""Set up the pdb2cif package."""
import setuptools


setuptools.setup(
    name="pdb2cif",
    version="1.1.0",
    description=(
        "This code reads and writes the old Protein Data Bank format."
    ),
    long_description=(
        "This code reads and writes the old PDB format "
        "(https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html) "
        "and allows conversion to/from mmCIF format (https://mmcif.wwpdb.org/)."
    ),
    python_requires=">=3.8",
    license="CC0-1.0",
    author="Nathan Baker",
    author_email="nathanandrewbaker@gmail.com",
    url="https://github.com/Electrostatics/pdb2cif",
    packages=setuptools.find_packages(),
    package_data={"": ["tests/data/*"]},
    install_requires=["requests", "pandas", "mmcif_pdbx"],
    extras_require={
        "dev": ["check-manifest"],
        "test": [
            "black",
            "coverage",
            "flake8",
            "pandas >= 1.0",
            "pytest",
            "testfixtures",
        ],
    },
    tests_require=[
        "pandas >= 1.0",
        "pytest",
        "testfixtures",
    ],
    test_suite="tests",
    # entry_points={"console_scripts": ["mvalue=osmolytes.main:console"]},
    keywords="science chemistry biophysics biochemistry",
    classifiers=[
        "Development Status :: 3 - Alpha",
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
