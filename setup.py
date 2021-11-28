import setuptools

with open("README.md", "r", encoding="utf-8") as readme:
    long_description = readme.read()

setuptools.setup(
    name="azapy",
    version="0.0.1",
    author="Mircea Marinescu",
    author_email="mircea.marinescu@outlook.com",
    description="Financial portfolio optimization algorithms",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    project_urls={
        "Bug Tracker": "",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPLv3 License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "azapy"},
    packages=setuptools.find_packages(where="azapy"),
    python_requires=">=3.8",
)
