import setuptools
import versioneer

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rayoptics",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author="Michael J Hayford",
    author_email="mjhoptics@gmail.com",
    description="Tools for image forming optical design and analysis",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    license="BSD-3-Clause",
    url="https://github.com/mjhoptics/ray-optics",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    keywords=['geometric optics', 'ray tracing', 'image forming optics',
              'paraxial optics', 'optical design', 'lens design',
              'aberrations', 'opd', 'psf'],
    install_requires=[
        "opticalglass",
        "numpy>=1.15.0",
        "scipy>=1.1.0",
        "matplotlib>=2.2.3",
        "json_tricks>=3.12.1",
        "pandas>=0.23.4",
        "attrs>=18.1.0",
        "transforms3d>=0.3.1"
        ],
    extras_require={
        'QtGUI':  ["pyqt5"],
    },
    entry_points={
        'gui_scripts': [
            'rayoptics = rayoptics.qtgui.rayopticsapp:main [QtGUI]',
        ],
    },
    data_files=[
        ('codev', ['tla_mapping.csv', 'test/*.seq', 'test/*.roa']),
        ('util', ['cie-cmf.txt', '*.csv']),
        ],
)
