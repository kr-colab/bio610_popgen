book: $(wildcard book/*.ipynb) $(wildcard book/_*.yml)
	jupyter-book build book/
	rm -rf docs/*
	cp -R book/_build/html/* docs/

clean: 
	jupyter-book clean book/
