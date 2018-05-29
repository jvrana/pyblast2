PIP=pip3


init:
	$(PIP) install pipenv --upgrade
	pipenv install --dev --skip-lock


build:
#	pipenv lock --requirements > requirements.txt
#	pipenv lock --dev --requirements > requirements_dev.txt
	docker build . --tag pyblast:testing

blast:
	pipenv run python blast_bin/install_blast.py $(EMAIL) $(PLATFORM)


test:
	pipenv run py.test


ci:
	pipenv run py.test -n 8 --boxed


testwithcoverage:
	pipenv run py.test --cov=pyblast --cov-report term-missing


pylint:
	pipenv run pylint pyblast

