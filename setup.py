import setuptools
#from azapy import __version__

with open("README.md", "r", encoding="utf-8") as readme:
    long_description = readme.read()

setuptools.setup(
    name="azapy",
    version="1.2.3",
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
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "."},
    packages=setuptools.find_packages(),
    python_requires=">=3.11",
    install_requires=[
          'numpy',
          'pandas',
          'scipy',
          'plotly',
          'matplotlib',
          'requests',
          'ecos',
          'pandas_market_calendars',
          'cvxopt',
          'ta',
          'yfinance',
          'statsmodels'
      ],
)
