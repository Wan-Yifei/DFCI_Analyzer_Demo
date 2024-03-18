from setuptools import setup, find_namespace_packages, find_packages

# Read requirements from requirements.txt
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='dfci_analysis_toolkit_demo',
    entry_points={
        'console_scripts': [
            'dfci-analysis = dfci_analysis:main'
        ]
    },
    py_modules=['dfci_analysis'],
    packages=find_namespace_packages(),
    install_requires=requirements,
    author='Yifei Wan',
    author_email='yifei.wan@gmail.com',
    description='Toolkit to analyze FASTQ/FASTA files',
    url='https://github.com/Wan-Yifei/DFCI_Analyzer_Demo.git',
)
