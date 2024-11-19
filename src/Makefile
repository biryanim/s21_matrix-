CC=g++
CFLAGS=-std=c++17 -Wall -Wextra -Werror
LIBS =  -lgtest -lgcov
LINUX = -lrt -lpthread -lm -lsubunit 
LFLAGS = -fprofile-arcs -ftest-coverage 	
OS = $(shell uname -s)


all: s21_matrix_oop.a test


s21_matrix_oop.a:clean
	$(CC) $(CC_FLAGS) -c s21_matrix_oop.cc -o ./s21_matrix_oop.o
	ar rcs s21_matrix_oop.a s21_matrix_oop.o
	ranlib s21_matrix_oop.a

# mytest:
# 	$(CC) $(CC_FLAGS) s21_matrix_oop.cc test.cc -o test
# 	./test

test: s21_matrix_oop.a
ifeq ($(OS), Darwin)
	$(CC) $(CFLAGS) unit_tests.cc s21_matrix_oop.a -o unit_tests $(LIBS) 
else
	$(CC) unit_tests.cc s21_matrix_oop.a -o unit_tests $(LIBS) $(LINUX) 
endif
	./unit_tests


gcov_report:
ifeq ($(OS), Darwin)
	$(CC) $(CFLAGS) $(LFLAGS) s21_*.cc unit_tests.cc -o test.out $(LIBS) 
else
	$(CC) $(CFLAGS) $(LFLAGS) s21_*.cc unit_tests.cc -o test.out $(LIBS) $(LINUX)
endif
	./test.out
	mkdir report
	gcovr *.gcda
	gcovr -b
	gcovr
	gcovr --html-details -o ./report/report.html
# lcov -t "test" --no-external -o test.info -c -d .
# genhtml -o ./report test.info
# open ./report/index.html

valgrind: test
	valgrind --leak-check=full --track-origins=yes --show-leak-kinds=all -s ./unit_tests

style:
	cp ../materials/linters/.clang-format ./.clang-format
	clang-format -i ./*.cc ./*.h
	clang-format -n ./*.cc ./*.h
	rm -rf .clang-format

clean:
	rm -rf  *.o *.a
	rm -rf *.gcda
	rm -rf *.gcno
	rm -rf *.info
	rm -rf *.gcov
	rm -rf report
	rm -rf gcov_test
	rm -rf test.out