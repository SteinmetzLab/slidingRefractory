import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open('requirements.txt') as f:
    require = [x.strip() for x in f.readlines() if not x.startswith('git+')]

setuptools.setup(
    name="slidingRP",
    version="1.1.0",
    author="Noam Roth",
    description="Quality metric from spike trains",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/SteinmetzLab/slidingRefractory",
    project_urls={
        "Bug Tracker": "https://github.com/SteinmetzLab/slidingRefractory/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=require,
    package_dir={'': 'python'},
    packages=setuptools.find_packages(where="python"),
#    package_data={'viewephys': ['raster.ui', 'nav_file.ui', 'viewephys.svg']},
    python_requires=">=3.8",
)
