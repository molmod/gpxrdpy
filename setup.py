import setuptools
setuptools.setup(
    name="gpxrdpy",
    version="0.0.1",
    author="Sander Borgmans",
    author_email="sander.borgmans@ugent.be",
    description="PXRD pattern calculator based on pyobjcryst",
    url="https://github.com/SanderBorgmans/gpxrdpy",
    packages=setuptools.find_packages(),
    scripts=['scripts/gpxrd.py', 'scripts/gpxrd_average.sh'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=['numpy>=1.5', 'scipy>=0.17.1', 'h5py>=2.0.0', 'matplotlib>1.0.0'],
    python_requires='>=3.6',
)
