import setuptools

with open("README.md", "r", encoding="utf-8") as readme:
    long_description = readme.read()

setuptools.setup(
    name="azapy",
    version="0.0.1",
    author="Mircea Marinescu",
    author_email="mircea.marinescu@outlook.com",
    description="Financial Portfolio Optimization Algorithms",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Mircea-MMXXI/azapy.git",
    project_urls={
        "Documentation": "https://azapy.readthedocs.io/en/latest",
        "Source": "https://github.com/Mircea-MMXXI/azapy",
        "Bug Tracker": "https://github.com/Mircea-MMXXI/azapy/issues",
    },
    license="GPLv3",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPLv3 License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "azapy"},
    packages=setuptools.find_packages(where="azapy"),
    python_requires=">=3.8",
)
