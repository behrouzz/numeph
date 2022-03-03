import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="numeph",
    version="0.3.0",
    author="Behrouz Safari",
    author_email="behrouz.safari@gmail.com",
    description="Convert JPL SPK ephemeris to numpy array",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/behrouzz/numeph",
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=["numeph"],
    include_package_data=True,
    install_requires=["numpy", "jplephem"],
    python_requires='>=3.4',
)
