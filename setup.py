import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="metapop",
    version="0.0.23",
    author="Kenji Gerhardt",
    author_email="kenji.gerhardt@gmail.com",
    description="A metagenomic statistic pipeline.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
	py_modules=["metapop"],
	include_package_data=True,
    python_requires='>=3',
	entry_points={
        "console_scripts": [
            "metapop=metapop.metapop_main:main",
        ]
    }
)