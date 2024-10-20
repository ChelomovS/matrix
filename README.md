# Matrix
Проектирование матрицы и поиск определителя

## Установка
Склонируйте репозиторий, перейдите в папку с репозиторием.

Сборка проекта:
```sh
cmake -B build
cmake --build build
```

## Использование 
Перейти в папку build:
```sh
cd build 
```

Для запуска программы:
```sh
./matrix
```
Для запуcка unit тестирования:
```sh
cd tests
./unit_tests
```

## End to end тестирование
Перейти в папку с тестами:
```sh
cd tests
```

Запустить end to end тестирование:
```sh
cd end_to_end
./end_to_end_test.sh
```