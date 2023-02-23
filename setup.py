import setuptools

setuptools.setup(name='priori',
      version='1.0.3',
      description="""Priori predicts transcription factor activity from RNA sequencing data using prior, literature-supported regulatory relationship information.""",
      author='Joseph Estabrook, William Yashar, & Emek Demir',
      author_email='yashar@ohsu.edu',
      entry_points={'console_scripts':['priori=priori.priori:main']},
      packages=setuptools.find_packages(where="priori"),
      package_dir={"": "priori"},
      include_package_data=True,
      url = 'https://github.com/ohsu-comp-bio/regulon-enrichment',
      classifiers = [
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: GLP-3.0",
         "Operating System :: OS Independent",
      ],
      python_requires = '==3.6.15',
      install_requires=[
          'numpy==1.19.5',
          'pandas==1.1.5',
          'scikit-learn==0.23.2',
          'scipy==1.5.3',
          'tqdm==4.64.0',
          'dill==0.3.5.1',
          'statsmodels==0.12.2'
        ]
     )
