from setuptools import setup

setup(
    name='mian',
    version='0.1.0',
    description="",
    long_description=(open('README.md').read()),
    url='https://github.com/tbj128/mian',
    install_requires=['Flask', 'Werkzeug', 'flask_login', 'Flask-Mail', 'scipy', 'pandas', 'h5py', 'rpy2', 'biom-format', 'scikit-bio', 'scikit-learn'],
    license='MIT',
    author='Tom Jin',
    author_email='',
    packages=['mian'],
    entry_points={
        'console_scripts': [
            'mian = gene.main:main',
        ],
    },
    include_package_data=True,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    zip_safe=False,
)
