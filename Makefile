PIP=pip3
.PHONY: docs  # necessary so it doesn't look for 'docs/makefile html'

init:
	pip install pip -U
	curl -sSL https://raw.githubusercontent.com/sdispater/poetry/master/get-poetry.py | python
	poetry self:update
	poetry install
	poetry run pre-commit install


clean:
	rm -rf dist
	rm -rf pip-wheel-metadata
	rm -rf docs
	rm -rf .pytest_cache


blast:
	poetry run python blast_bin/install_blast.py $(EMAIL) $(PLATFORM)


testwithcoverage:
	poetry run py.test --cov=pyblast --cov-report term-missing


test:
	poetry run python -m pytest


lint:
	poetry run pylint -E pydent


docs:
	@echo "Updating docs"

	# copy README.md to README.rst format for Sphinx documentation
	# we can comment this out if we do not want to include the README.md in the sphinx documentation

	rm -rf docs
	cd docsrc && poetry run make html
	find docs -type f -exec chmod 444 {} \;
	@echo "\033[95m\n\nBuild successful! View the docs homepage at docs/html/index.html.\n\033[0m"

	touch docs/.nojekyll
	open ./docs/index.html


format:
	poetry run keats version up
	poetry run keats changelog up
	poetry run black pyblast tests


ci:
	poetry run py.test -n 8 --boxed


lock:
	poetry run upver
	poetry update



build:
	poetry run upver
	poetry build


release:
	sh scripts/release.sh


klocs:
	find . -name '*.py' | xargs wc -l
