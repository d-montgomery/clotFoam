.PHONY: help docker-build docker-test docker-clean docker-shell

help:
	@echo "clotFoam - OpenFOAM 12 Docker Testing"
	@echo ""
	@echo "Available targets:"
	@echo "  make docker-build  - Build Docker image with OpenFOAM 12"
	@echo "  make docker-test   - Run compilation and testing in Docker"
	@echo "  make docker-shell  - Open interactive shell in Docker container"
	@echo "  make docker-clean  - Remove Docker images and containers"

docker-build:
	@echo "Building Docker image..."
	docker build -t clotfoam:openfoam12 -f docker/Dockerfile .

docker-test: docker-build
	@echo "Running tests in Docker..."
	docker run --rm \
		clotfoam:openfoam12 \
		/bin/bash /home/ofuser/docker/run-test.sh

docker-shell: docker-build
	@echo "Opening interactive shell..."
	docker run --rm -it \
		clotfoam:openfoam12 \
		/bin/bash -c "source /usr/lib/openfoam/openfoam2412/etc/bashrc && /bin/bash"

docker-clean:
	@echo "Cleaning Docker artifacts..."
	docker rmi clotfoam:openfoam12 || true
	docker system prune -f
