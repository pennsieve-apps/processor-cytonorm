.PHONY: help test fixture run image clean

SERVICE_NAME ?= processor-cytonorm
IMAGE ?= $(SERVICE_NAME):dev
DEMO ?= ../workflow-immune-cytometry/data/zenodo-15388817
EVENTS ?= 20000
PYTHON ?= python3
export DOCKER_DEFAULT_PLATFORM = linux/amd64

.DEFAULT: help

help:
	@echo "Make Help for $(SERVICE_NAME)"
	@echo ""
	@echo "make test     - unit tests (tests/run_tests.R) on a synthetic fixture; uses a local Rscript"
	@echo "                with the packages installed, else runs them in the image"
	@echo "make fixture  - build data/fixture/{fcs,metadata} from the demo data (EVENTS per file, 0 = all;"
	@echo "                needs numpy + flowio: pip install -r requirements-dev.txt)"
	@echo "make run      - build the image and run it with docker compose over data/fixture into data/output"
	@echo "make clean    - remove data/output"

test:
	@if command -v Rscript >/dev/null 2>&1 && Rscript -e 'quit(status = !requireNamespace("flowCore", quietly = TRUE))' >/dev/null 2>&1; then \
		Rscript tests/run_tests.R; \
	else \
		docker build -q -t $(IMAGE) . >/dev/null && \
		docker run --rm -v "$(CURDIR):/src:ro" -w /src --entrypoint Rscript $(IMAGE) tests/run_tests.R; \
	fi

fixture:
	$(PYTHON) tests/make_fixture.py --source $(DEMO) --out data/fixture --events $(EVENTS)

image:
	docker build -t $(IMAGE) .

run: clean
	@test -d data/fixture/fcs || { echo "no data/fixture: run 'make fixture' first"; exit 1; }
	docker compose -f docker-compose.yml down --remove-orphans
	docker compose -f docker-compose.yml build
	docker compose -f docker-compose.yml up --exit-code-from processor
	docker compose -f docker-compose.yml down --remove-orphans
	@find data/output -type f | sort

clean:
	rm -rf data/output
