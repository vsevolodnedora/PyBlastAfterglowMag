import setuptools

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setuptools.setup(
    name='PyBlastAfterglowMag',
    author='Vsevolod Nedora',
    author_email='vsevolod.nedora@gmail.com',
    description='Python package that accompanies C++ PyBlastAfterglowMag code',
    keywords='afterglow, GRB',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/PartumSomnia/PyBlastAfterglowMag',
    project_urls={
        'Documentation': 'https://github.com/PartumSomnia/PyBlastAfterglowMag',
        'Bug Reports':
            'https://github.com/PartumSomnia/PyBlastAfterglowMag/issues',
        'Source Code': 'https://github.com/PartumSomnia/PyBlastAfterglowMag',
        # 'Funding': '',
        # 'Say Thanks!': '',
    },
    package_dir={'': 'src'},
    packages=setuptools.find_packages(where='src'),
    classifiers=[
        # see https://pypi.org/classifiers/
        'Development Status :: 5 - Production/Stable',

        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    # install_requires=['Pillow'],
    extras_require={
        'dev': ['check-manifest'],
        # 'test': ['coverage'],
    },
    # entry_points={
    #     'console_scripts': [  # This can provide executable scripts
    #         'run=examplepy:main',
    # You can execute `run` in bash to run `main()` in src/examplepy/__init__.py
    #     ],
    # },
)