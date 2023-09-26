from setuptools import setup

setup(
    name="mtz",
    version="0.1",
    description="Simple mtz-reading library",
    author="Nicholas Devenish",
    author_email="ndevenish@gmail.com",
    package_dir={"": "src"},
    packages=["mtz"],
    entry_points={
        "console_scripts": [
            "mtzshow = mtz.mtz:run",
        ]
    },
)
