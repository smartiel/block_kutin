from setuptools import setup
from setuptools_rust import Binding, RustExtension

setup(
    name="block-kutin",
    version="0.1.0",
    author="Simon Martiel and Timothée Goubault de Brugière",
    rust_extensions=[RustExtension("block_kutin.block_kutin", binding=Binding.PyO3)],
    packages=["block_kutin"],
    # rust extensions are not zip safe, just like C-extensions.
    zip_safe=False,
)
