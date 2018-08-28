import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rayoptics",
    version="0.1.5.dev1",
    author="Michael J Hayford",
    author_email="mjhoptics@gmail.com",
    description="Tools for image forming optical design and analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="BSD-3-Clause",
    url="https://github.com/mjhoptics/ray-optics",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "json_tricks",
        "pandas",
        "pyqt5",
        "attr"
        ]
)
