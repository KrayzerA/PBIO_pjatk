import random
from collections import Counter

# Cel programu:
# Program generuje losową sekwencję DNA w formacie FASTA, zgodnie z wymaganiami użytkownika.
# Zawiera możliwość wstawienia imienia w losowe miejsce w sekwencji, a także oblicza statystyki zawartości nukleotydów.

# Program jest przydatny w bioinformatyce do generowania przykładowych sekwencji DNA, które można wykorzystać do testów lub edukacji.

# Funkcja generująca losową sekwencję DNA
def generate_random_sequence(length):
    nucleotides = ['A', 'C', 'G', 'T']  # Lista możliwych nukleotydów
    return ''.join(random.choice(nucleotides) for _ in range(length))  # Generowanie losowej sekwencji


# Funkcja do wstawiania imienia w losowe miejsce
def insert_name_in_sequence(sequence, name):
    insert_position = random.randint(0, len(sequence) - 1)  # Losowe miejsce w sekwencji na imię
    # ORIGINAL (created by AI)
    # sequence = sequence[:insert_position] + name + sequence[insert_position + len(name):]
    # MODIFIED (fixed: nie poprawnie wstawianie imienia do sekwencji, byly tracone nukleotydy przy wstawianu imienia w srodek)
    return sequence[:insert_position] + name.lower() + sequence[insert_position:]  # Wstawiamy imię do sekwencji


# Funkcja zapisująca sekwencję do pliku FASTA
def save_to_fasta(name, description, sequence):
    file_name = f"{name}.fasta"  # Nazwa pliku
    with open(file_name, "w") as fasta_file:
        fasta_file.write(f">{name} {description}\n")  # Nagłówek pliku FASTA
        fasta_file.write(sequence + "\n")  # Zapisanie sekwencji DNA do pliku
    print(f"Zapisano sekwencję do pliku {file_name}")  # Potwierdzenie zapisu


# Funkcja do obliczania statystyk
def calculate_statistics(sequence):
    counts = Counter(sequence)  # Zliczamy wystąpienia każdego nukleotydu
    total_length = len(sequence)  # Całkowita długość sekwencji

    # Obliczamy procentową zawartość każdego nukleotydu
    percentages = {nucleotide: (counts[nucleotide] / total_length) * 100 for nucleotide in 'ACGT'}

    cg_percentage = percentages['C'] + percentages['G'] # procentowy udział C + G

    return percentages, cg_percentage


def generate_dna_sequence(length, name, description, my_name):
    # ORIGINAL (created by AI)
    # nucleotides = ['A', 'C', 'G', 'T']
    # sequence = ''.join(random.choice(nucleotides) for _ in range(length))
    # MODIFIED (refactored)
    sequence = generate_random_sequence(length)

    # ORIGINAL (created by AI)
    # insert_position = random.randint(0, length - 1)
    # sequence = sequence[:insert_position] + my_name + sequence[insert_position + len(my_name):]
    # MODIFIED (refactored)
    sequence_with_name = insert_name_in_sequence(sequence, my_name)

    # ORIGINAL (created by AI)
    # count_a = sequence.count('A')
    # count_c = sequence.count('C')
    # count_g = sequence.count('G')
    # count_t = sequence.count('T')
    # total_length = len(sequence)
    # percentage_a = (count_a / total_length) * 100
    # percentage_c = (count_c / total_length) * 100
    # percentage_g = (count_g / total_length) * 100
    # percentage_t = (count_t / total_length) * 100
    # cg_percentage = ((count_c + count_g) / total_length) * 100
    # MODIFIED (refactored and optimised)
    percentages, cg_percentage = calculate_statistics(sequence)

    # ORIGINAL (created by AI)
    # print(f"Statystyki sekwencji:")
    # print(f"A: {percentage_a:.1f}%")
    # print(f"C: {percentage_c:.1f}%")
    # print(f"G: {percentage_g:.1f}%")
    # print(f"T: {percentage_t:.1f}%")
    # print(f"%CG: {cg_percentage:.1f}")
    # MODIFIED (refactored)
    # Wyświetlanie statystyk
    print("Statystyki sekwencji:")
    for nucleotide, percentage in percentages.items():
        print(f"{nucleotide}: {percentage:.1f}%")
    print(f"%CG: {cg_percentage:.1f}")

    # ORIGINAL (created by AI)
    # file_name = f"{name}.fasta"
    # with open(file_name, "w") as fasta_file:
    #     fasta_file.write(f">{name} {description}\n")
    #     fasta_file.write(sequenceWithName + "\n")
    # MODIFIED (refactored)
    save_to_fasta(name, description, sequence_with_name)


if __name__ == "__main__":
    try:
        length = int(input("Podaj długość sekwencji: "))  # Długość sekwencji
        name = input("Podaj ID sekwencji: ")  # ID sekwencji
        description = input("Podaj opis sekwencji: ")  # Opis sekwencji
        my_name = input("Podaj imię: ")  # Imię do wstawienia w sekwencji

        # Wywołanie funkcji do generowania sekwencji
        generate_dna_sequence(length, name, description, my_name)

    except ValueError as e:
        print(f"Błąd: {e}")
