from collections import defaultdict
import matplotlib.pyplot as plt


def load_quality_profile(filename: str) -> dict[int, dict[str, dict]]:
    data = defaultdict(lambda: defaultdict(dict))
    with open(filename, "r") as f:
        lines = f.readlines()

    for i in range(0, len(lines), 2):
        if i + 1 >= len(lines):
            break

        score_line = lines[i].strip().split("\t")
        freq_line = lines[i + 1].strip().split("\t")

        if len(score_line) < 2 or len(freq_line) < 2:
            continue

        try:
            base = score_line[0]
            position = int(score_line[1])
            scores = list(map(int, score_line[2:]))
            frequencies = list(map(int, freq_line[2:]))
        except (ValueError, IndexError):
            continue

        cumulative = []
        running_total = 0

        for freq in frequencies:
            running_total += freq
            cumulative.append(running_total)

        try:
            total = cumulative[-1]
        except IndexError:
            total = 0

        data[position][base] = {
            "scores": scores,
            "cumulative": cumulative,
            "total": total
        }

    return {pos: dict(bases) for pos, bases in data.items()}


if __name__ == "__main__":

    profile = load_quality_profile(r"D:\Code\DNA-Sequence-Simulator\src\simulators\illumina\profiles\ecoli_hiseqR1.txt")

    position_scores = {}

    for idx in profile:
        total_weight = 0
        weighted_sum = 0

        for base in profile[idx]:
            scores = profile[idx][base]["scores"]
            cumulative = profile[idx][base]["cumulative"]
            total = profile[idx][base]["total"]

            if not cumulative or total == 0:
                continue

            weights = [cumulative[0]] + [
                cumulative[i] - cumulative[i - 1]
                for i in range(1, len(cumulative))
            ]

            weighted_sum += sum(s * w for s, w in zip(scores, weights))
            total_weight += total

        if total_weight > 0:
            position_scores[idx] = weighted_sum / total_weight

    x = sorted(position_scores)
    y = [position_scores[i] for i in x]

    plt.ylim(36, 40)
    plt.plot(x, y)
    plt.xlabel("Read Position")
    plt.ylabel("Mean Quality Score")
    plt.title("Mean Quality Score by Position")
    plt.tight_layout()
    plt.show()