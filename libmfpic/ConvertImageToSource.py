import numpy as np
from PIL import Image

def grayscale_intensity_2d(image_filepath, desired_width):
  image = Image.open(image_filepath)
  original_width, original_height = image.size
  desired_height = original_height * desired_width // original_width
  image = image.resize((desired_width, desired_height), resample = Image.Resampling.LANCZOS)
  new_width, new_height = image.size
  gray = image.convert("L")
  intensity = np.array(gray)
  normed_intensity = intensity / 255.0
  return normed_intensity, new_width, new_height

def writeSourceFile(image_filepath, average_number_density, desired_num_cells_wide, cell_length):
  image_intensity, num_cells_wide, num_cells_high = grayscale_intensity_2d(image_filepath, desired_num_cells_wide)
  average_darkness = np.mean(image_intensity)
  number_density_intensity_scalar = average_number_density / average_darkness

  num_cells = num_cells_high * num_cells_wide
  domain_width = num_cells_wide * cell_length
  domain_height = num_cells_high * cell_length

  source_code = f"""
if (x[0] <= {cell_length} and x[1] >= {domain_height - cell_length}) {{
  number_density = {number_density_intensity_scalar * image_intensity[0, 0]};
}}
  """
  for (y_pixel_from_top, x_pixel), pixel_intensity in np.ndenumerate(image_intensity):
    source_code += f"""
else if (x[0] <= {(x_pixel + 1) * cell_length} and x[1] >= {domain_height - (y_pixel_from_top + 1) * cell_length}) {{
  number_density = {number_density_intensity_scalar * pixel_intensity};
}}
    """

  with open("source.txt", 'w') as source_file:
    source_file.write(source_code)

  print(f"I am {num_cells_wide} cells wide and {num_cells_high} high")

if __name__ == "__main__":
  writeSourceFile("mfpic.png", 1.0e16, 113, 1.0e-6)
