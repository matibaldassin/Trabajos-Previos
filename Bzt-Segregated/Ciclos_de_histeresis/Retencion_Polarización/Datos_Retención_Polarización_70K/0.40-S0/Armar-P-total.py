pasos = iter(range(10000, 1500001, 10000))

with open("P-total.txt", "w") as out:
    for i in range(1, 16):
        with open(f"P-total-files/P-total-{i}", "r") as f:
            data = f.read().splitlines()
            if i != 1:
                data = data[2:]
            for a in range(len(data)):
                if i == 1 and a < 2:
                    out.write(f"\t{data[a]}\n")
                else:
                    out.write(f"{next(pasos)}\t{'  '.join(data[a].split()[1:])}\n")
