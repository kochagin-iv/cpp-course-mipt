Написать шаблонный класс `Finite<int N>`, реализующий концепцию конечного поля из N элементов. Должны поддерживаться арифметические операции, кроме деления, а в случае простого N и деление тоже. Элемент поля Finite<N> должно быть можно сконструировать от int путем взятия остатка от деления на N (в математическом смысле).

После этого, используя ранее написанный класс рациональных чисел, написать класс `Matrix` с тремя шаблонными параметрами: unsigned M, unsigned N, typename Field=Rational. (По умолчанию берётся поле рациональных чисел, но можно создать матрицу и над конечным полем.)

Матрицы должны поддерживать следующие операции:

Проверка на равенство: операторы `==` и `!=`.<br>
Конструктор по умолчанию, создающий единичную матрицу. Для неквадратных матриц конструктор по умолчанию не требуется.<br>
Конструктор, создающий матрицу из vector<vector<T> >. Должно быть можно создать матрицу из vector<vector<int> >.<br>
Сложение, вычитание, операторы `+=`, `-=`. Сложение и вычитание матриц несоответствующих размеров не должно компилироваться.<br>
Умножение на число.<br>
Умножение матриц, работающие за max(M,N,K)**3. Для квадратных матриц должен поддерживаться еще и оператор *=. Попытка перемножить матрицы несоответствующих размеров должна приводить к ошибке компиляции.<br>
Метод `det()`, возвращающий определитель матрицы за O(N**3). Взятие определителя от неквадратной матрицы не должно компилироваться.<br>
Метод `transposed()`, возвращающий транспонированную матрицу.<br>
Метод `rank()` - вычислить ранг матрицы.<br>
Метод `trace()` - вычислить след матрицы.<br>
Методы `inverted()` и `invert()` - вернуть обратную матрицу и обратить данную матрицу.<br>
Методы `getRow(unsigned)` и `getColumn(unsigned)`, возвращающие `std::vector<Field>` из соответствующих значений.<br>
К матрице должен быть дважды применим оператор `[]`, причём это должно работать как для неконстантных, так и для константных матриц. В первом случае содержимое матрицы должно быть можно таким способом поменять.<br>
Другие способы изменения содержимого матрицы, кроме описанных выше, должны отсутствовать. Однако не запрещается реализовать дополнительные методы для выполнения каких-либо иных алгебраических операций или для удобства работы, если по названию и сигнатуре этих методов будет без комментариев понятно их действие.<br>
Квадратные матрицы размера N должно быть можно объявлять всего с одним обязательным шаблонным параметром: SquareMatrix<N>.<br>
