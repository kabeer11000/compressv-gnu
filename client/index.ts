import fs from "fs";
import sharp from "sharp";
import axios from "axios";
import path from "path";

async function main() {
  const inputPath = path.resolve("input.jpg"); // change this to your image path
  const targetWidth = 64;
  const targetHeight = 64;

  try {
    // Resize and convert image to raw grayscale
    const rawBuffer = await sharp(inputPath)
      .resize(targetWidth, targetHeight)
      .grayscale()
      .raw()
      .toBuffer();

    if (rawBuffer.length !== 4096) {
      console.error(`Expected 4096 bytes, got ${rawBuffer.length}`);
      process.exit(1);
    }

    console.log("Sending image data to /compress_image...");

    const response = await axios.post(
      "http://localhost:18080/compress_image",
      rawBuffer,
      {
        headers: {
          "Content-Type": "application/octet-stream",
          "Content-Length": rawBuffer.length,
        },
        responseType: "arraybuffer", // or "json" depending on your backend
      }
    );

    console.log("Server responded with status:", response.status);
    // Optional: save response if it's binary
    fs.writeFileSync("output_response.bin", Buffer.from(response.data));
    console.log("Saved response to output_response.bin");

  } catch (err) {
    console.error("Error:", err);
  }
}

main();
