"""Exact finite verification for the accompanying research paper.

The proofs in ``paper/manuscript.tex`` do not depend on this program.  This
script exhausts every 3 by 3 array over {empty, 1, 2, 3}, retains the monotone
partial matrices, and checks the finite instances listed in Appendix A.  It
also checks small cases of the two explicit constructions.

Only the Python standard library is used.  All non-integral calculations use
``fractions.Fraction``; there is no floating-point certification.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
from typing import DefaultDict, Iterable, Sequence


Cell = tuple[int, int]
Triple = tuple[int, int, int]
Matrix = list[list[int]]  # 0 denotes an empty cell.


def choose2(value: int) -> int:
    """Return binomial(value, 2)."""

    return value * (value - 1) // 2


def is_monotone(matrix: Sequence[Sequence[int]]) -> bool:
    """Test the three defining monotonicity conditions."""

    rows = len(matrix)
    columns = len(matrix[0])

    for row in matrix:
        values = [value for value in row if value]
        if any(left >= right for left, right in zip(values, values[1:])):
            return False

    for column in range(columns):
        values = [matrix[row][column] for row in range(rows) if matrix[row][column]]
        if any(top >= bottom for top, bottom in zip(values, values[1:])):
            return False

    by_label: DefaultDict[int, list[Cell]] = defaultdict(list)
    for row in range(rows):
        for column in range(columns):
            label = matrix[row][column]
            if label:
                by_label[label].append((row, column))

    for cells in by_label.values():
        cells.sort()
        if any(c >= next_c for (_, c), (_, next_c) in zip(cells, cells[1:])):
            return False
    return True


def position_maps(
    matrix: Sequence[Sequence[int]],
) -> tuple[list[dict[int, int]], list[dict[int, int]], dict[int, list[Cell]]]:
    """Return row positions, column positions, and cells grouped by label."""

    rows = len(matrix)
    columns = len(matrix[0])
    row_positions = [dict() for _ in range(rows)]
    column_positions = [dict() for _ in range(columns)]
    by_label: DefaultDict[int, list[Cell]] = defaultdict(list)
    for row in range(rows):
        for column in range(columns):
            label = matrix[row][column]
            if label:
                row_positions[row][label] = column
                column_positions[column][label] = row
                by_label[label].append((row, column))
    return row_positions, column_positions, dict(by_label)


def forcing_data(
    matrix: Sequence[Sequence[int]],
) -> tuple[dict[Cell, tuple[int, ...]], dict[Cell, str], dict[tuple[Cell, Cell], int], int]:
    """Build the forcing graph and assert its normal-form properties."""

    rows = len(matrix)
    columns = len(matrix[0])
    row_positions, column_positions, by_label = position_maps(matrix)

    forcing: dict[Cell, tuple[int, ...]] = {}
    orientation: dict[Cell, str] = {}
    for row in range(rows):
        for column in range(columns):
            if matrix[row][column]:
                continue
            labels = tuple(sorted(set(row_positions[row]) & set(column_positions[column])))
            forcing[row, column] = labels
            signs: set[str] = set()
            for label in labels:
                is_left = row_positions[row][label] < column
                is_below = column_positions[column][label] > row
                assert is_left == is_below
                signs.add("L" if is_left else "R")
            assert len(signs) <= 1
            if signs:
                orientation[row, column] = signs.pop()

    edges: dict[tuple[Cell, Cell], int] = {}
    adjacency: DefaultDict[Cell, list[tuple[Cell, int]]] = defaultdict(list)
    for label, cells in by_label.items():
        cells.sort()
        for first, second in combinations(cells, 2):
            row, column = first
            next_row, next_column = second
            left = (row, next_column)
            right = (next_row, column)
            assert orientation[left] == "L"
            assert orientation[right] == "R"
            edge = (left, right)
            assert edge not in edges
            edges[edge] = label
            adjacency[left].append((right, label))
            adjacency[right].append((left, label))

    # Degrees and ordered neighbourhoods.
    for hole, labels in forcing.items():
        if not labels:
            continue
        incident = sorted(adjacency[hole])
        assert len(incident) == len(labels)
        assert sorted(label for _, label in incident) == list(labels)
        for ((r, c), label), ((next_r, next_c), next_label) in zip(
            incident, incident[1:]
        ):
            assert r < next_r and c < next_c and label < next_label

    # Every colour class is an order-preserving matching.
    colour_edges: DefaultDict[int, list[tuple[Cell, Cell]]] = defaultdict(list)
    for edge, label in edges.items():
        colour_edges[label].append(edge)
    for same_colour in colour_edges.values():
        assert len({left for left, _ in same_colour}) == len(same_colour)
        assert len({right for _, right in same_colour}) == len(same_colour)
        for (left, right), (next_left, next_right) in combinations(same_colour, 2):
            if left[0] < next_left[0] and left[1] < next_left[1]:
                assert right[0] < next_right[0] and right[1] < next_right[1]
            if next_left[0] < left[0] and next_left[1] < left[1]:
                assert next_right[0] < right[0] and next_right[1] < right[1]

    # Every four-cycle is rainbow.
    left_vertices = [hole for hole in adjacency if orientation[hole] == "L"]
    for left, next_left in combinations(left_vertices, 2):
        neighbours = dict(adjacency[left])
        next_neighbours = dict(adjacency[next_left])
        common = sorted(set(neighbours) & set(next_neighbours))
        for right, next_right in combinations(common, 2):
            colours = {
                neighbours[right],
                neighbours[next_right],
                next_neighbours[right],
                next_neighbours[next_right],
            }
            assert len(colours) == 4

    expected_edges = sum(choose2(len(cells)) for cells in by_label.values())
    assert len(edges) == expected_edges
    assert sum(len(labels) for labels in forcing.values()) == 2 * expected_edges
    energy = sum(choose2(len(labels)) for labels in forcing.values())
    assert energy == sum(choose2(len(adjacency[hole])) for hole in adjacency)
    return forcing, orientation, edges, energy


def inversion_cells(matrix: Sequence[Sequence[int]]) -> set[frozenset[Cell]]:
    """Enumerate strict inversions as unordered pairs of physical cells."""

    filled = [
        (row, column, matrix[row][column])
        for row in range(len(matrix))
        for column in range(len(matrix[0]))
        if matrix[row][column]
    ]
    inversions: set[frozenset[Cell]] = set()
    for (row, column, label), (next_row, next_column, next_label) in combinations(
        filled, 2
    ):
        if row < next_row and column < next_column and label > next_label:
            inversions.add(frozenset(((row, column), (next_row, next_column))))
        if next_row < row and next_column < column and next_label > label:
            inversions.add(frozenset(((row, column), (next_row, next_column))))
    return inversions


def saturated_images(
    matrix: Sequence[Sequence[int]],
    forcing: dict[Cell, tuple[int, ...]],
    orientation: dict[Cell, str],
) -> set[frozenset[Cell]]:
    """Construct both sides of the saturated-inversion bijection and compare."""

    row_positions, column_positions, _ = position_maps(matrix)
    images: set[frozenset[Cell]] = set()
    event_count = 0
    for (row, column), labels in forcing.items():
        for smaller, larger in combinations(labels, 2):
            event_count += 1
            if orientation[row, column] == "L":
                first = (row, row_positions[row][larger])
                second = (column_positions[column][smaller], column)
            else:
                first = (column_positions[column][larger], column)
                second = (row, row_positions[row][smaller])
            images.add(frozenset((first, second)))
    assert len(images) == event_count

    enumerated: set[frozenset[Cell]] = set()
    for inversion in inversion_cells(matrix):
        first, second = sorted(inversion)
        row, column = first
        next_row, next_column = second
        assert row < next_row and column < next_column
        labels = {matrix[row][column], matrix[next_row][next_column]}
        saturated_crosses = [
            cross
            for cross in ((row, next_column), (next_row, column))
            if set(forcing.get(cross, ())) >= labels
        ]
        assert len(saturated_crosses) <= 1
        if saturated_crosses:
            enumerated.add(inversion)
    assert images == enumerated
    return images


def cyclic_representations(
    matrix: Sequence[Sequence[int]],
) -> list[tuple[Matrix, dict[Cell, Triple]]]:
    """Return the three coordinate-as-label representations and cell maps."""

    n = len(matrix)
    triples = [
        (row, column, matrix[row][column] - 1)
        for row in range(n)
        for column in range(n)
        if matrix[row][column]
    ]
    answer: list[tuple[Matrix, dict[Cell, Triple]]] = []
    for label_axis in range(3):
        free_axes = [axis for axis in range(3) if axis != label_axis]
        representation = [[0] * n for _ in range(n)]
        cell_to_triple: dict[Cell, Triple] = {}
        for triple in triples:
            cell = (triple[free_axes[0]], triple[free_axes[1]])
            representation[cell[0]][cell[1]] = triple[label_axis] + 1
            cell_to_triple[cell] = triple
        assert is_monotone(representation)
        answer.append((representation, cell_to_triple))
    return answer


def aligned_pair_count(matrix: Sequence[Sequence[int]]) -> int:
    """Count pairs for which all three strict coordinate signs agree."""

    triples = [
        (row, column, matrix[row][column] - 1)
        for row in range(len(matrix))
        for column in range(len(matrix[0]))
        if matrix[row][column]
    ]
    count = 0
    for first, second in combinations(triples, 2):
        signs = [(a > b) - (a < b) for a, b in zip(first, second)]
        if 0 not in signs and (all(sign == 1 for sign in signs) or all(sign == -1 for sign in signs)):
            count += 1
    return count


def check_exact_slack(
    original: Sequence[Sequence[int]],
    representations: Sequence[Sequence[Sequence[int]]],
    forcing_by_axis: Sequence[dict[Cell, tuple[int, ...]]],
    energies: Sequence[int],
) -> None:
    """Check the pair partition and exact slack identity over the rationals."""

    n = len(original)
    m = sum(bool(value) for row in original for value in row)
    holes = n * n - m
    if holes == 0:
        return

    first_counts: list[int] = []
    forcing_variances: list[Fraction] = []
    inversion_counts: list[int] = []
    for matrix, forcing in zip(representations, forcing_by_axis):
        _, _, by_label = position_maps(matrix)
        first_count = sum(len(cells) * (len(cells) - 1) for cells in by_label.values())
        assert first_count == sum(len(labels) for labels in forcing.values())
        mean = Fraction(first_count, holes)
        variance = sum(
            (Fraction(len(labels)) - mean) ** 2 for labels in forcing.values()
        )
        first_counts.append(first_count)
        forcing_variances.append(variance)
        inversion_counts.append(len(inversion_cells(matrix)))

    tied_pairs = Fraction(sum(first_counts), 2)
    aligned_pairs = aligned_pair_count(original)
    assert Fraction(choose2(m)) == tied_pairs + aligned_pairs + sum(inversion_counts)

    unsaturated = sum(
        inversions - energy
        for inversions, energy in zip(inversion_counts, energies)
    )
    assert unsaturated >= 0
    right_side = (
        sum(Fraction(first_count * first_count, 2 * holes) for first_count in first_counts)
        + sum(forcing_variances, Fraction()) / 2
        + aligned_pairs
        + unsaturated
    )
    assert Fraction(choose2(m)) == right_side

    for first_count, matrix in zip(first_counts, representations):
        _, _, by_label = position_maps(matrix)
        degrees = [len(by_label.get(label, ())) for label in range(1, n + 1)]
        coordinate_variance = sum((Fraction(degree) - Fraction(m, n)) ** 2 for degree in degrees)
        assert Fraction(first_count) == Fraction(m * m, n) - m + coordinate_variance


def check_endpoint_exclusion(matrix: Sequence[Sequence[int]], alphabet_size: int) -> None:
    """Check endpoint exclusion for every row subset of size at least two."""

    rows = len(matrix)
    columns = len(matrix[0])
    row_degrees = [sum(bool(value) for value in row) for row in matrix]
    for size in range(2, rows + 1):
        for row_set in combinations(range(rows), size):
            common_columns = {
                column
                for column in range(columns)
                if all(matrix[row][column] for row in row_set)
            }
            label_sets = [{value for value in matrix[row] if value} for row in row_set]
            common_labels = set.intersection(*label_sets)
            assert size * len(common_labels) + len(common_columns) <= columns
            assert len(common_labels) + size * len(common_columns) <= alphabet_size
            assert len(common_labels) + len(common_columns) <= min(
                row_degrees[row] for row in row_set
            )


def pressure_records(
    matrix: Sequence[Sequence[int]],
) -> tuple[dict[Cell, tuple[tuple[int, ...], set[Cell]]], list[dict[int, int]], list[dict[int, int]]]:
    """Build pressure grids and check the triangular-hole estimate."""

    n = len(matrix)
    row_positions, column_positions, _ = position_maps(matrix)
    records: dict[Cell, tuple[tuple[int, ...], set[Cell]]] = {}
    for row in range(n):
        for column in range(n):
            if matrix[row][column]:
                continue
            labels = tuple(sorted(set(row_positions[row]) & set(column_positions[column])))
            if not labels:
                continue
            arm_rows = tuple(column_positions[column][label] for label in labels)
            arm_columns = tuple(row_positions[row][label] for label in labels)
            assert (
                all(arm_row > row for arm_row in arm_rows)
                and all(arm_column < column for arm_column in arm_columns)
            ) or (
                all(arm_row < row for arm_row in arm_rows)
                and all(arm_column > column for arm_column in arm_columns)
            )
            grid = {(arm_row, arm_column) for arm_row in arm_rows for arm_column in arm_columns}
            hole_count = sum(matrix[r][c] == 0 for r, c in grid)
            assert 2 * hole_count >= len(labels) * (len(labels) + 1)
            records[row, column] = labels, grid
    return records, row_positions, column_positions


def is_two_comparable(triples: Iterable[Triple]) -> bool:
    """Test pairwise two-comparability."""

    triples = list(triples)
    for first, second in combinations(triples, 2):
        smaller = sum(a < b for a, b in zip(first, second))
        larger = sum(a > b for a, b in zip(first, second))
        if max(smaller, larger) < 2:
            return False
    return True


def check_recursive_pressure(matrix: Sequence[Sequence[int]]) -> None:
    """Check both encodings and all projection bounds for every physical cell."""

    records, row_positions, column_positions = pressure_records(matrix)
    n = len(matrix)
    m = sum(bool(value) for row in matrix for value in row)
    for physical_row in range(n):
        for physical_column in range(n):
            signed_centres: list[list[tuple[int, int, tuple[int, ...]]]] = [[], []]
            for (row, column), (labels, grid) in records.items():
                if (physical_row, physical_column) not in grid:
                    continue
                if row < physical_row and physical_column < column:
                    signed_centres[0].append((row, column, labels))
                elif physical_row < row and column < physical_column:
                    signed_centres[1].append((row, column, labels))
                else:
                    raise AssertionError("pressure-grid quadrant failure")

            total_centres = 0
            total_load = 0
            for centres in signed_centres:
                first_encoding: list[Triple] = []
                second_encoding: list[Triple] = []
                incidences: set[Triple] = set()
                selected_first_cells: set[Cell] = set()
                selected_second_cells: set[Cell] = set()
                for row, column, labels in centres:
                    a_label = matrix[physical_row][column]
                    b_label = matrix[row][physical_column]
                    assert a_label in labels and b_label in labels
                    p_value = row_positions[row][a_label]
                    q_value = column_positions[column][b_label]
                    first_encoding.append((row, column, p_value))
                    second_encoding.append((row, column, q_value))
                    selected_first_cells.add((row, p_value))
                    selected_second_cells.add((q_value, column))
                    incidences.update((row, column, label) for label in labels)

                assert is_two_comparable(first_encoding)
                assert is_two_comparable(second_encoding)
                assert len(selected_first_cells) == len(centres)
                assert len(selected_second_cells) == len(centres)
                assert len(centres) <= m

                projection_12 = {(row, column) for row, column, _ in incidences}
                projection_13 = {(row, label) for row, _, label in incidences}
                projection_23 = {(column, label) for _, column, label in incidences}
                assert len(projection_12) <= m
                assert len(projection_13) <= m
                assert len(projection_23) <= m
                assert len(incidences) ** 2 <= (
                    len(projection_12) * len(projection_13) * len(projection_23)
                )
                assert len(incidences) ** 2 <= m**3
                total_centres += len(centres)
                total_load += len(incidences)

            assert total_centres <= 2 * m
            assert total_load * total_load <= 4 * m**3


def biclique_matrix(s: int, t: int) -> Matrix:
    """Return the matrix from Proposition 3.5."""

    n = s * t
    matrix = [[0] * n for _ in range(n)]
    for a in range(s):
        for b in range(t):
            label = a + s * b + 1
            matrix[a][b] = label
            matrix[s + b][t + a] = label
    assert is_monotone(matrix)
    return matrix


def sharp_pressure_matrix(s: int, t: int) -> tuple[Matrix, Cell, int]:
    """Return the matrix and claimed load from Proposition 7.5."""

    block = s + t
    n = s + s * block + 1
    matrix = [[0] * n for _ in range(n)]
    first_pressure_column = 1 + s * block
    bottom_row = n - 1

    for i in range(s):
        matrix[i][0] = i + 1
        for u in range(block):
            matrix[i][1 + i * block + u] = s + u + 1

    for j in range(s):
        column = first_pressure_column + j
        for u in range(block):
            matrix[s + j * block + u][column] = u + 1
        matrix[bottom_row][column] = s + t + j + 1

    assert is_monotone(matrix)
    return matrix, (bottom_row, 0), s * s * (t + 2)


def check_explicit_constructions() -> None:
    """Check the biclique and sharp-pressure constructions in small cases."""

    for s, t in ((2, 2), (2, 3), (3, 3), (3, 4)):
        matrix = biclique_matrix(s, t)
        forcing, _, edges, _ = forcing_data(matrix)
        left = [(a, t + a) for a in range(s)]
        right = [(s + b, b) for b in range(t)]
        assert all((x, y) in edges for x in left for y in right)
        assert all(len(forcing[x]) == t for x in left)

    for s, t in ((1, 1), (2, 2), (3, 3), (4, 6)):
        matrix, cell, expected_load = sharp_pressure_matrix(s, t)
        records, _, _ = pressure_records(matrix)
        actual_load = sum(
            len(labels) for labels, grid in records.values() if cell in grid
        )
        assert actual_load == expected_load
        assert matrix[cell[0]][cell[1]] == 0
        check_recursive_pressure(matrix)


def main() -> None:
    """Run the exhaustive search and all construction checks."""

    arrays_examined = 0
    monotone_matrices = 0
    cyclic_representation_checks = 0
    total_energy = 0

    for flat in product(range(4), repeat=9):
        arrays_examined += 1
        matrix = [list(flat[3 * row : 3 * row + 3]) for row in range(3)]
        if not is_monotone(matrix):
            continue
        monotone_matrices += 1
        check_recursive_pressure(matrix)

        physical_saturated_images: list[set[frozenset[Triple]]] = []
        representations = cyclic_representations(matrix)
        forcing_by_axis: list[dict[Cell, tuple[int, ...]]] = []
        energies: list[int] = []
        representation_matrices: list[Matrix] = []
        for representation, cell_to_triple in representations:
            cyclic_representation_checks += 1
            forcing, orientation, _, energy = forcing_data(representation)
            check_endpoint_exclusion(representation, 3)
            images = saturated_images(representation, forcing, orientation)
            physical_saturated_images.append(
                {frozenset(cell_to_triple[cell] for cell in image) for image in images}
            )
            forcing_by_axis.append(forcing)
            energies.append(energy)
            representation_matrices.append(representation)
            total_energy += energy

        assert all(
            first.isdisjoint(second)
            for first, second in combinations(physical_saturated_images, 2)
        )
        m = sum(bool(value) for row in matrix for value in row)
        assert sum(energies) <= choose2(m)
        check_exact_slack(matrix, representation_matrices, forcing_by_axis, energies)

    assert arrays_examined == 4**9
    assert monotone_matrices == 712
    assert cyclic_representation_checks == 3 * monotone_matrices
    check_explicit_constructions()

    print(f"arrays examined: {arrays_examined}")
    print(f"monotone order-three matrices: {monotone_matrices}")
    print(f"cyclic representations checked: {cyclic_representation_checks}")
    print(f"aggregate second-energy count checked: {total_energy}")
    print("explicit biclique cases checked: (2,2), (2,3), (3,3), (3,4)")
    print("explicit pressure cases checked: (1,1), (2,2), (3,3), (4,6)")
    print("all finite verification assertions passed")


if __name__ == "__main__":
    main()
