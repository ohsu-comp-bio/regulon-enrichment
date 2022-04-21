import setuptools

setuptools.setup(name='regulon-enrichment',
      version='0.4.a0',
      description="""This package leverages pathway information and gene expression data to produce
        regulon-based protein activity scores""",
      long_description = """regulon-enrichment leverages pathway information and gene expression data to produce
        regulon-based protein activity scores""",
      author='Joseph Estabrook',
      author_email='estabroj@ohsu.edu',
      entry_points={'console_scripts':['enrich=enricher.enrich:main']},
      packages=setuptools.find_packages(
          exclude=["enricher.tests.*", "enricher.tests", "enricher.plot.*","enricher.plot"]),
      data_files=[('enricher',
                   ['data/PathwayCommons9.All.hgnc.sif.gz', 'data/secondary_intx_regulon.pkl']
                   )],
      package_data={'enricher':['data/PathwayCommons9.All.hgnc.sif.gz', 'data/secondary_intx_regulon.pkl']},
      url = 'https://github.com/JEstabrook/regulon-enrichment',
      classifiers = [
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
      ],
      python_requires = '>=3.6',
      install_requires=[
          'numpy>=1.17.3',
          'pandas>=0.25.3',
          'scikit-learn>=0.21.3',
          'scipy>=1.3.1',
          'tqdm>=4.38.0',
          'dill>=0.3.1.1'
        ]
     )
