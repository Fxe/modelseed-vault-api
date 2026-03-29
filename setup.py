from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="modelseed_vault",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "requests>=2.25.0",  # For making HTTP requests
        "pydantic>=1.8.0",   # For data validation and settings management
        "typing-extensions>=3.7.4",  # For additional type hints
        "python-dotenv>=0.19.0",  # For environment variable management
        "aiohttp>=3.8.0",  # For async HTTP requests
        "biopython>=1.79",  # For sequence handlings
        "modelseedpy>=0.1.0",  # For ModelSEED API
        "cobra>=0.20.0",  # For COBRA APIs
        "pymongo",
        "neo4j",
        "lxml",
    ],
    author="Filipe Liu",
    author_email="fliu@anl.gov",
    description="ModelSEED Vault Python API",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/fxe/modelseed_annotation_api",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Natural Language :: English",
    ],
    python_requires=">=3.7",
)
