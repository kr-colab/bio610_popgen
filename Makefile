book: book/*
	jupyter-book build book/
	rm -rf docs/*
	cp -R book/_build/html/* docs/

clean: 
	jupyter-book clean book/
