from setuptools import setup, find_packages


setup(
    name="entrezclient",
    version="0.1.0",
    description="Entrez client",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Pieter Provoost",
    author_email="pieterprovoost@gmail.com",
    url="https://github.com/pieterprovoost/entrez-client",
    packages=find_packages(),
    python_requires=">=3.6"
)