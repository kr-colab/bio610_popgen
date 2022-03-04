.PHONY : book

book : docs/.book_build

docs/.book_build : $(wildcard book/*.ipynb) $(wildcard book/_*.yml)
	jupyter-book build book/
	rm -rf docs/*
	cp -R book/_build/html/* docs/
	touch $@

clean : 
	jupyter-book clean book/
	rm docs/.book_build
