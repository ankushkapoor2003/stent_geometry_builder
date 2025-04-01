from setuptools import setup, find_packages

setup(
    name='sld_interface',
    version='0.0.2',
    description='Integration of Solidworks API with Matlab using pywin32',
    url="https://github.com/ankushkapoor2003/stent_geometry_builder",
    author='Ankush',
    author_email='ankushkapoor2003au@gmail.com',
    license='All Rights Reserved',
    python_requires='==3.9.18',
    packages=find_packages(),
    install_requires=[
        'pywin32==306',
        'numpy==1.26.3',
    ],
    zip_safe=False
)