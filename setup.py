from setuptools import setup

pkg = 'reliontomotools'

exec(open(f'{pkg}/_version.py').read())

setup(
      name=pkg,
      version=__version__,
      packages=[pkg],
      entry_points={
          "console_scripts": [
              f"warptomo2relion = {pkg}.warptomo2relion:warpTomo2RelionProgram"
              ]
          },
      )
