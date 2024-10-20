import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def generate_pairs(alignments):
    # Сначала создадим список для хранения всех пар
    pairs = []

    # Определим количество строк в выравнивании
    num_sequences = len(alignments)

    # Пройдем по каждой позиции (столбцу) в выравнивании
    for i in range(len(alignments[0])):
        # Получим нуклеотиды в текущей позиции
        column = [alignments[j][i] for j in range(num_sequences)]

        # Генерируем все возможные пары
        for x in column:
            for y in column:
                pairs.append((x, y))

    return pairs


def count_pairs(alignments):
    pair_counts = {}

    # Получаем количество строк и столбцов
    num_rows = len(alignments)
    num_cols = len(alignments[0])

    # Проходим по каждому столбцу
    for col in range(num_cols):
        # Собираем нуклеотиды в текущем столбце
        nucleotides = []
        for row in range(num_rows):
            nucleotide = alignments[row][col]
            if nucleotide != '-':
                nucleotides.append(nucleotide)

        # Подсчитываем пары
        for i in range(len(nucleotides)):
            for j in range(i, len(nucleotides)):
                n1 = nucleotides[i]
                n2 = nucleotides[j]

                # Создаем кортеж для пары
                pair = (n1, n2) if n1 <= n2 else (n2, n1)

                # Увеличиваем счетчик для пары
                if pair in pair_counts:
                    pair_counts[pair] += 1
                else:
                    pair_counts[pair] = 1

    return pair_counts


def calculate_frequencies(alignments):
    # Инициализация счетчика для нуклеотидов
    nucleotide_counts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
    total_nucleotides = 0

    # Проходим по каждому столбцу выравнивания
    for col in range(len(alignments[0])):
        for row in alignments:
            nucleotide = row[col]
            if nucleotide in nucleotide_counts:  # Проверяем, что это нуклеотид
                nucleotide_counts[nucleotide] += 1
                total_nucleotides += 1

    # Вычисляем частоты
    frequencies = {nuc: count / total_nucleotides for nuc, count in nucleotide_counts.items() if total_nucleotides > 0}

    return frequencies


# Примеры использования функций

alignments = [
    "AGCTACGTGTCGCTGAATCTATGACT",
    "-GCTA-GAGCA-AGGCAACTGCATCT",
    "A-CTG-CACCC-ATGAACCTCGCGCT",
    "A-CTG-CACCC-ATGAACCTCTCGCT",
    "A-CTG-CACCC-ATGAACCTCTCGCT",
    "A-CTG-CACCC-ATGAACCTCTCACT",
    "A-CTG-CACCC-ATGAACCTCTCACT"
]

# Генерация пар
pairs_result = generate_pairs(alignments)
print(f"Количество пар: {len(pairs_result)}")

# Подсчет пар
pair_counts = count_pairs(alignments)
print("Подсчет пар:")
print(pair_counts)

# Вычисление частот нуклеотидов
freqs = calculate_frequencies(alignments)
print("Частоты:")
for nuc, freq in freqs.items():
    print(f"{nuc}: {freq:.4f}")


def align(top_seq, bottom_seq, gap_penalty, blosum_matrix):
    # Инициализация размеров матрицы
    n = len(top_seq)
    m = len(bottom_seq)

    # Создание матрицы A размером (n+1) x (m+1)
    A = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

    # Заполнение первой строки и первого столбца
    for i in range(n + 1):
        A[i][0] = -i * gap_penalty
    for j in range(m + 1):
        A[0][j] = -j * gap_penalty

    # Заполнение матрицы A
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = blosum_matrix[top_seq[i - 1]][bottom_seq[j - 1]]
            match = A[i - 1][j - 1] + match_score
            delete = A[i - 1][j] - gap_penalty
            insert = A[i][j - 1] - gap_penalty

            A[i][j] = max(match, delete, insert)

    return A


def get_alignment(top_seq, bottom_seq, A, gap_penalty, blosum_matrix):
    aligned_top = []
    aligned_bottom = []

    i, j = len(top_seq), len(bottom_seq)

    while i > 0 and j > 0:
        score_current = A[i][j]
        score_diagonal = A[i - 1][j - 1]
        score_up = A[i][j - 1]
        score_left = A[i - 1][j]

        if score_current == score_diagonal + blosum_matrix[top_seq[i - 1]][bottom_seq[j - 1]]:
            aligned_top.append(top_seq[i - 1])
            aligned_bottom.append(bottom_seq[j - 1])
            i -= 1
            j -= 1
        elif score_current == score_left - gap_penalty:
            aligned_top.append(top_seq[i - 1])
            aligned_bottom.append('-')
            i -= 1
        else:  # score_current == score_up - gap_penalty
            aligned_top.append('-')
            aligned_bottom.append(bottom_seq[j - 1])
            j -= 1

    while i > 0:
        aligned_top.append(top_seq[i - 1])
        aligned_bottom.append('-')
        i -= 1

    while j > 0:
        aligned_top.append('-')
        aligned_bottom.append(bottom_seq[j - 1])
        j -= 1

    aligned_top.reverse()
    aligned_bottom.reverse()

    return ''.join(aligned_top), ''.join(aligned_bottom)

result_matrix = align(top_seq, bottom_seq, gap_penalty, blosum_matrix)

# Вывод матрицы
print("Matrix:")
for row in result_matrix:
    print(row)

# Получение и вывод выравненных последовательностей
aligned_sequences = get_alignment(top_seq, bottom_seq, result_matrix, gap_penalty, blosum_matrix)
print("nAligned Sequences:")
print(aligned_sequences[0])
print(aligned_sequences[1])


def align(top_seq, bottom_seq, gap_penalty, blosum_matrix):
    n = len(top_seq) + 1
    m = len(bottom_seq) + 1

    # Инициализация матрицы
    A = [[0] * m for _ in range(n)]

    # Заполнение первой строки и первого столбца
    for i in range(n):
        A[i][0] = -i * gap_penalty
    for j in range(m):
        A[0][j] = -j * gap_penalty

    # Заполнение матрицы
    for i in range(1, n):
        for j in range(1, m):
            match_score = blosum_matrix[top_seq[i - 1]][bottom_seq[j - 1]]
            A[i][j] = max(
                A[i - 1][j - 1] + match_score,  # совпадение/несоответствие
                A[i - 1][j] - gap_penalty,  # пропуск в нижней последовательности
                A[i][j - 1] - gap_penalty  # пропуск в верхней последовательности
            )

    return A


def get_alignment(top_seq, bottom_seq, A, gap_penalty, blosum_matrix):
    aligned_top = []
    aligned_bottom = []

    i, j = len(top_seq), len(bottom_seq)

    while i > 0 and j > 0:
        score_current = A[i][j]
        score_diagonal = A[i - 1][j - 1]
        score_up = A[i][j - 1]
        score_left = A[i - 1][j]

        # Проверяем, был ли это совпадение или несоответствие
        if score_current == score_diagonal + blosum_matrix[top_seq[i - 1]][bottom_seq[j - 1]]:
            aligned_top.append(top_seq[i - 1])
            aligned_bottom.append(bottom_seq[j - 1])
            i -= 1
            j -= 1
        elif score_current == score_left - gap_penalty:
            aligned_top.append(top_seq[i - 1])
            aligned_bottom.append('-')
            i -= 1
        else:  # score_current == score_up - gap_penalty
            aligned_top.append('-')
            aligned_bottom.append(bottom_seq[j - 1])
            j -= 1

    # Если остались незаконченные последовательности
    while i > 0:
        aligned_top.append(top_seq[i - 1])
        aligned_bottom.append('-')
        i -= 1

    while j > 0:
        aligned_top.append('-')
        aligned_bottom.append(bottom_seq[j - 1])
        j -= 1

    # Реверсируем строки, так как мы добавляли символы в обратном порядке
    aligned_top.reverse()
    aligned_bottom.reverse()

    return ''.join(aligned_top), ''.join(aligned_bottom)

sm = align(top_seq, bottom_seq, gap_penalty, blosum_matrix)
aligned_sequences = get_alignment(top_seq, bottom_seq, sm, gap_penalty, blosum_matrix)
print("Aligned sequences for first example:")
print(aligned_sequences[0])
print(aligned_sequences[1])

def create_blosum_matrix(scores, nucleotides):
    """Создает симметричную матрицу BLOSUM из заданных оценок."""
    matrix = {nuc: {nuc2: 0 for nuc2 in nucleotides} for nuc in nucleotides}
    for (nuc1, nuc2), score in scores.items():
        matrix[nuc1][nuc2] = score
        matrix[nuc2][nuc1] = score  # Симметричная матрица
    return matrix


def visualize_blosum_matrix(matrix, nucleotides):
    """Визуализирует матрицу BLOSUM с помощью тепловой карты."""
    data = np.array([[matrix[x][y] for y in nucleotides] for x in nucleotides])
    plt.figure(figsize=(10, 8))
    sns.heatmap(data, xticklabels=nucleotides, yticklabels=nucleotides, annot=True, cmap="coolwarm")
    plt.title("Матрица BLOSUM")
    plt.show()


def init(n, m, penalty):
    """Инициализирует матрицу для алгоритма выравнивания."""
    c = [[0] * (m + 1) for _ in range(n + 1)]

    # Заполнение первой строки
    for j in range(1, m + 1):
        c[0][j] = c[0][j - 1] - penalty

    # Заполнение первого столбца
    for i in range(1, n + 1):
        c[i][0] = c[i - 1][0] - penalty

    return c


def fill_matrix(a, b, blosum, penalty):
    """Заполняет матрицу выравнивания на основе последовательностей и матрицы BLOSUM."""
    n = len(a)
    m = len(b)

    # Инициализация матрицы
    c = init(n, m, penalty)

    # Заполнение остальной части матрицы
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = c[i - 1][j - 1] + blosum[a[i - 1]][b[j - 1]]
            delete = c[i - 1][j] - penalty
            insert = c[i][j - 1] - penalty
            c[i][j] = max(match, delete, insert)

    return c


def get_new_score(up, left, middle, s_score, gap_penalty):
    """Вычисляет новое значение для ячейки матрицы."""
    match = middle + s_score  # Совпадение
    delete = up - gap_penalty  # Удаление
    insert = left - gap_penalty  # Вставка
    return max(match, delete, insert)

# Создание и визуализация матрицы BLOSUM
blosum_matrix = create_blosum_matrix(scores, nucleotides)
visualize_blosum_matrix(blosum_matrix, nucleotides)

matrix = fill_matrix(a, b, blosum_matrix, penalty)
for row in matrix:
    print(row)

new_score = get_new_score(up_value, left_value, middle_value, s_score, gap_penalty)
print("Новое значение:", new_score)  # Вывод нового значения для матрицы


