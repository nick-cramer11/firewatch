import os
import re
import zipfile

def unzip_and_rename_all(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for filename in os.listdir(input_dir):
        match = re.match(r'era5_(\d{2})_(\d{4})\.nc$', filename)
        if not match:
            continue  # skip files that don't match pattern

        month, year = match.groups()
        zip_path = os.path.join(input_dir, filename)

        try:
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                for file in zip_ref.namelist():
                    orig_name = os.path.basename(file)
                    name, ext = os.path.splitext(orig_name)
                    new_name = f"{name}_{month}_{year}{ext}"
                    output_path = os.path.join(output_dir, new_name)

                    with zip_ref.open(file) as source, open(output_path, "wb") as target:
                        target.write(source.read())

                print(f"✅ Unzipped and renamed contents of {filename}")
        except zipfile.BadZipFile:
            print(f"⚠️ Skipping {filename} — not a valid ZIP archive.")

# Example usage:
unzip_and_rename_all("data", "data/unzipped")
