def remove_values_by_index(filename, indices_to_remove):
    """
    Reads a file containing space-separated values, removes the values at specified indices (except the first line),
    and writes the modified content back to the file.

    :param filename: The name of the file to process.
    :param indices_to_remove: A list of indices (0-based) to remove from each line (except the first line).
    """
    # Read the original file content
    with open(filename, 'r') as f:
        lines = f.readlines()

    # The first line will not be modified
    modified_lines = [lines[17]]  # Keep the first line as is

    # Process each subsequent line by removing the values at the specified indices
    for line in lines[18:]:
        values = line.split()  # Split the line into values
        # Remove the values at specified indices (make sure to sort indices in reverse order for safe removal)
        modified_values = [value for idx, value in enumerate(values) if idx not in indices_to_remove]
        modified_lines.append(" ".join(modified_values) + "\n")

    # Write the modified content back to the file
    with open(filename, 'w') as f:
        f.writelines(modified_lines)

    print(f"Successfully modified {filename} by removing values at indices: {indices_to_remove}, except the first line.")

# Example usage:
filename = "p.txt"  # Replace with your actual file
indices_to_remove = [1,2,3]  # Remove values at these indices (0-based, except for the first line)

remove_values_by_index(filename, indices_to_remove)

