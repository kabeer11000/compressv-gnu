import fs from "fs";
import sharp from "sharp";

async function main() {
  const rawPath = "output_response.bin"; // JSON file with raw_hex field
  const outputImagePath = "output.png";

  const rawBuffer = fs.readFileSync(rawPath);
  const json = JSON.parse(rawBuffer.toString());

  if (!json.raw_hex) {
    console.error("JSON does not contain 'raw_hex' field");
    return;
  }

  // Decode hex string to Buffer
  const imgBuffer = Buffer.from(json.raw_hex, "hex");

  if (imgBuffer.length !== 4096) {
    console.error(`Expected 4096 bytes for image data, got ${imgBuffer.length}`);
    return;
  }

  await sharp(imgBuffer, {
    raw: {
      width: 64,
      height: 64,
      channels: 1, // grayscale
    },
  })
    .png()
    .toFile(outputImagePath);

  console.log(`Saved as ${outputImagePath}`);
}

main();
