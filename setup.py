from distutils.core import setup
setup(name='genomeToTranscriptMapper',
      version='0.1',
      install_requires=[
                      'numpy',                     
                      ],
      packages=['genomeToTranscriptMapper'],
      entry_points={
        'console_scripts': [
            'genomeToTranscriptMapper=genomeToTranscriptMapper.genomeToTranscriptMapper:run'
        ]
      }
      )