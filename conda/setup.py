from setuptools import setup

# @todo install a python file that imports sps from the conda lib folder...

setup(
    name='sps',
    version='0.2.0',
    description='O(1) region count queries using sparse prefix sums',
    url='https://github.com/MarkusRainerSchmidt/libSps',
    author='Markus Schmidt',
    author_email='markus.rainer.schmidt@gmail.com',
    license='MIT',
    packages=["sps"],
    package_dir={"": "pip_install"},
    package_data={"sps": ["libsps.cpython-36m-x86_64-linux-gnu.so"]},
    #ext_modules=[Extension('sps', libraries = ['sps'], sources=["../src/tree.cpp"])],
    install_requires=[
        # @todo currently done via conda environment
    ],

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
)