import setuptools
#from azapy import __version__

with open("README.md", "r", encoding="utf-8") as readme:
    long_description = readme.read()

setuptools.setup(
    name="azapy",
    version="1.2.4",
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
          'numpy>=1.26.0',
          'pandas>=2.2.0',
          'scipy>=1.12.0',
          'plotly>=5.19.0',
          'matplotlib>=3.8.0',
          'requests>=2.31.0',
          'ecos>=2.0.0',
          'pandas_market_calendars>=4.4.0',
          'cvxopt>=1.3.2',
          'ta>=0.11.0',
          'yfinance>=0.2.37',
          'statsmodels>=0.14.0'
      ],
)
