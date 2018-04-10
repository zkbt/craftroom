# craftroom
A good craftroom is strewn with pipe cleaners, bits of yarn, scraps of fabric, colorful paper, and at least one hot glue gun. When it comes time to make a puppet or a toy or a Halloween costume, the ingredients you need are probably sitting on a shelf somewhere. This is a collection of useful code snippets and tools that hopes to serve the same messy purpose for various astronomy projects.

You should be able to install it either by cloning the repository and running `python setup.py install` from its main directory, or directly by running `pip install git+https://github.com/zkbt/craftroom.git`.

If you want to be able to modify the code yourself, please also feel free to fork/clone this repository onto your own computer and install directly from that editable package. For example, this might look like:
```
git clone https://github.com/zkbt/craftroom.git
cd craftroom/
pip install -e .
```
This will link the installed version of the `craftroom` package to your local repository. Changes you make to the code in the repository should be reflected in the version Python sees when it tries to `import craftroom`.
