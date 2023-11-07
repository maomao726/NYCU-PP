CXX=g++
CXXFLAGS = -O3 -pthread -ffast-math

APP_NAME=pi

${APP_NAME}: ${APP_NAME}.cpp
	${CXX} $< ${CXXFLAGS} -o $@.out

clean:
	rm *.out
