import setuptools

with open("README.md", "r", encoding="utf-8") as readme:
    long_description = readme.read()

setuptools.setup(
    name="azapy",
    version="1.2.6",
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
    python_requires=">=3.14",
    install_requires=[
          'numpy>=2.4.4',
          'pandas>=3.0.2',
          'scipy>=1.17.1',
          'plotly>=6.6.0',
          'matplotlib>=3.10.8',
          'requests>=2.31.1',
          'ecos>=2.0.14',
          'pandas_market_calendars>=5.3.2',
          'cvxopt>=1.3.3',
          'ta>=0.11.0',
          'yfinance>=1.3.0',
          'statsmodels>=0.14.6',
          'jinja2>=3.1.6',
          'nbformat>=5.10.4'
      ],
)
