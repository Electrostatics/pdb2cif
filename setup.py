"""Set up the pdb2cif package."""
import setuptools


setuptools.setup(
    name="pdb_old_format",
    version="0.0.1",
    description=(
        "This code reads and writes the old Protein Data Bank format."
    ),
    long_description=(
        "This code reads and writes the old PDB format "
        "(https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html)."
    ),
    python_requires=">=3.8",
    license="CC0-1.0",
    author="Nathan Baker",
    author_email="nathanandrewbaker@gmail.com",
    url="https://github.com/Electrostatics/pdb2cif",
    packages=setuptools.find_packages(),
    package_data={"": ["tests/data/*"]},
    install_requires=["requests"],
    tests_require=["pytest", "pandas"],
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
