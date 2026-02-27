import csv
import json
from datetime import datetime


WHO_VARIANTS = {
    "Alpha",
    "Beta",
    "Gamma",
    "Delta",
    "Epsilon",
    "Eta",
    "Iota",
    "Kappa",
    "Lambda",
    "Mu",
    "Omicron",
}


def parse_date(s: str):
    return datetime.strptime(s, "%Y-%m-%d")


def canonical_variant(name: str) -> str:
    start = name.find("(")
    end = name.find(")")
    if start != -1 and end != -1 and end > start + 1:
        inner = name[start + 1 : end]
        inner = inner.split(",")[0].strip()
        return inner or "Other"
    return name


def to_date_str(dt: datetime) -> str:
    return dt.strftime("%Y-%m-%d")


def main():
    workspace = "/Users/gcohn/Documents/data-vis-a3"
    cases_path = f"{workspace}/COVID_US_cases.csv"
    per_country_path = f"{workspace}/perCountryData.json"
    output_path = f"{workspace}/vis1_interpolated_daily_and_monthly_percentages.json"

    # Daily cases: drop last row, clamp negatives to 0.
    with open(cases_path, "r", encoding="utf-8") as f:
        rows_raw = list(csv.DictReader(f))
    rows = rows_raw[:-1] if rows_raw else []
    daily_cases = []
    for d in rows:
        date_str = d.get("date", "")
        if not date_str:
            continue
        date = parse_date(date_str)
        new_confirmed = max(0.0, float(d.get("new_confirmed", 0) or 0))
        daily_cases.append({"date": date, "dateStr": date_str, "newConfirmed": new_confirmed})
    daily_cases.sort(key=lambda x: x["date"])
    if not daily_cases:
        raise RuntimeError("No daily cases found in COVID_US_cases.csv")

    with open(per_country_path, "r", encoding="utf-8") as f:
        per_country = json.load(f)

    # Use only first USA distribution.
    usa_timepoints = []
    canonical_set = set()
    usa_distribution_used = False

    regions = per_country.get("regions", []) if isinstance(per_country, dict) else []
    for region in regions:
        if usa_distribution_used:
            break
        for entry in region.get("distributions", []):
            if usa_distribution_used:
                break
            if entry.get("country") != "USA" or not isinstance(entry.get("distribution"), list):
                continue
            usa_distribution_used = True
            for tp in entry["distribution"]:
                week = tp.get("week")
                total_seq = float(tp.get("total_sequences", 0) or 0)
                if not week or total_seq <= 0:
                    continue
                date = parse_date(week)
                counts = tp.get("cluster_counts", {}) or {}
                canonical_counts = {}
                for raw_name, raw_count in counts.items():
                    count = float(raw_count or 0)
                    if count <= 0:
                        continue
                    c = canonical_variant(str(raw_name))
                    c_lower = c.lower()
                    who_match = None
                    for w in WHO_VARIANTS:
                        if w.lower() == c_lower:
                            who_match = w
                            break
                    key = who_match if who_match else "non_who"
                    canonical_counts[key] = canonical_counts.get(key, 0.0) + count

                sum_canonical = sum(canonical_counts.values())
                if sum_canonical < total_seq:
                    canonical_counts["unknown"] = canonical_counts.get("unknown", 0.0) + (total_seq - sum_canonical)

                perc_by_variant = {}
                for c, val in canonical_counts.items():
                    if val <= 0:
                        continue
                    perc_by_variant[c] = val / total_seq
                    canonical_set.add(c)

                if perc_by_variant:
                    usa_timepoints.append({"date": date, "percByVariant": perc_by_variant})

    usa_timepoints.sort(key=lambda x: x["date"])
    canonical_set.add("unknown")
    canonical_set.add("non_who")

    # Interpolate to daily percentages.
    daily_perc = {}
    if not usa_timepoints:
        for dc in daily_cases:
            daily_perc[dc["date"]] = {"unknown": 1.0}
    else:
        first_tp_date = usa_timepoints[0]["date"]
        last_tp_date = usa_timepoints[-1]["date"]
        idx = 0

        for dc in daily_cases:
            t = dc["date"]
            perc_map = {}

            if t < first_tp_date:
                perc_map = {"unknown": 1.0}
            elif t > last_tp_date:
                perc_map = dict(usa_timepoints[-1]["percByVariant"])
            else:
                while idx + 1 < len(usa_timepoints) and t > usa_timepoints[idx + 1]["date"]:
                    idx += 1
                a = usa_timepoints[idx]
                b = usa_timepoints[min(idx + 1, len(usa_timepoints) - 1)]
                if a["date"] == b["date"]:
                    perc_map = dict(a["percByVariant"])
                else:
                    span = (b["date"] - a["date"]).total_seconds()
                    alpha = (t - a["date"]).total_seconds() / span if span else 0.0
                    keys = set(a["percByVariant"].keys()) | set(b["percByVariant"].keys())
                    for k in keys:
                        p0 = a["percByVariant"].get(k, 0.0)
                        p1 = b["percByVariant"].get(k, 0.0)
                        perc_map[k] = p0 + (p1 - p0) * alpha

                s = sum(perc_map.values())
                if s < 0.999:
                    perc_map["unknown"] = perc_map.get("unknown", 0.0) + (1.0 - s)

            daily_perc[t] = perc_map

    # Aggregate monthly variant case counts.
    monthly_map = {}
    for dc in daily_cases:
        cases = dc["newConfirmed"]
        if cases <= 0:
            continue
        month_key = dc["dateStr"][:7]
        bucket = monthly_map.get(month_key)
        if bucket is None:
            bucket = {"total": 0.0, "byVariant": {}}
            monthly_map[month_key] = bucket
        bucket["total"] += cases
        perc_map = daily_perc.get(dc["date"], {"unknown": 1.0})
        for v, frac in perc_map.items():
            if not frac or frac <= 0:
                continue
            add = cases * frac
            if add <= 0:
                continue
            bucket["byVariant"][v] = bucket["byVariant"].get(v, 0.0) + add

    # Output structures requested:
    # 1) interpolated daily variant values
    # 2) monthly percentage values
    interpolated_daily_variant_values = []
    for dc in daily_cases:
        perc_map = daily_perc.get(dc["date"], {})
        interpolated_daily_variant_values.append(
            {
                "date": dc["dateStr"],
                "new_confirmed": round(dc["newConfirmed"], 6),
                "variant_fractions": {k: round(v, 10) for k, v in sorted(perc_map.items()) if v > 0},
            }
        )

    monthly_percentage_values = []
    for month_key in sorted(monthly_map.keys()):
        bucket = monthly_map[month_key]
        total = bucket["total"]
        by_variant_cases = {k: v for k, v in bucket["byVariant"].items() if v > 0}
        by_variant_percentages = {}
        if total > 0:
            by_variant_percentages = {k: (v / total) for k, v in by_variant_cases.items()}
        monthly_percentage_values.append(
            {
                "month": month_key,
                "monthly_new_cases": round(total, 6),
                "variant_case_counts": {k: round(v, 6) for k, v in sorted(by_variant_cases.items())},
                "variant_percentages": {k: round(v, 10) for k, v in sorted(by_variant_percentages.items())},
            }
        )

    output = {
        "metadata": {
            "source_cases_csv": "COVID_US_cases.csv",
            "source_variants_json": "perCountryData.json",
            "usa_distribution_used": "first USA only",
            "who_variants": sorted(WHO_VARIANTS),
            "notes": [
                "Daily variant fractions are linearly interpolated between USA timepoints.",
                "Dates before first USA timepoint are set to unknown=1.",
                "Dates after last USA timepoint use the last USA timepoint percentages.",
                "Unknown receives remainder if interpolated fractions sum to < 0.999.",
            ],
        },
        "interpolated_daily_variant_values": interpolated_daily_variant_values,
        "monthly_percentage_values": monthly_percentage_values,
    }

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)

    print(f"Wrote {output_path}")
    print(f"daily rows: {len(interpolated_daily_variant_values)}")
    print(f"monthly rows: {len(monthly_percentage_values)}")


if __name__ == "__main__":
    main()
