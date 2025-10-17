using System.Numerics;
using System.Windows;
using System.Windows.Media;
using System.Windows.Media.Imaging;

namespace ColorGradient;

/// <summary>
/// Main window for the Color Gradient application.
/// Demonstrates proper color interpolation in linear RGB space and conversion to sRGB for display.
/// </summary>
public partial class MainWindow : Window
{
    /// <summary>
    /// Width of the rendered image in pixels.
    /// </summary>
    private const int Width = 512;
    
    /// <summary>
    /// Height of the rendered image in pixels.
    /// </summary>
    private const int Height = 256;
    
    /// <summary>
    /// Writeable bitmap used for rendering the gradient.
    /// </summary>
    private readonly WriteableBitmap _writeableBitmap;
    
    /// <summary>
    /// Initializes a new instance of the MainWindow class.
    /// Sets up the bitmap and renders the color gradient.
    /// </summary>
    public MainWindow()
    {
        InitializeComponent();
        _writeableBitmap = new WriteableBitmap(Width, Height, 96, 96, PixelFormats.Bgra32, null);
        Canvas.Source = _writeableBitmap;
        DrawGradientLinearRgbThenToSrgb();
    }
    
    /// <summary>
    /// Performs linear interpolation between two 3D vectors.
    /// This is used for interpolating color values in linear RGB space.
    /// </summary>
    /// <param name="a">The starting vector (at t=0).</param>
    /// <param name="b">The ending vector (at t=1).</param>
    /// <param name="t">The interpolation parameter, typically in range [0, 1].</param>
    /// <returns>The interpolated vector: a + (b - a) * t</returns>
    public static Vector3 Lerp(Vector3 a, Vector3 b, float t)
    {
        return a + (b - a) * t;
        
    }
    
    /// <summary>
    /// Converts a linear RGB color component to an 8-bit sRGB value.
    /// </summary>
    /// <param name="linear">Linear RGB component value (typically in range [0, 1]).</param>
    /// <returns>8-bit sRGB value in range [0, 255].</returns>
    /// <remarks>
    /// The sRGB transfer function applies a piecewise transformation:
    /// - For small values (≤ 0.0031308): linear scaling by 12.92
    /// - For larger values: power function with gamma ≈ 2.4
    /// This matches the standard sRGB specification (IEC 61966-2-1).
    /// </remarks>
    private static byte ToSrgb8(float linear)
    {
        if (linear <= 0f) return 0;
        if (linear >= 1f) return 255;

        float srgb = linear <= 0.0031308f
            ? 12.92f * linear
            : 1.055f * MathF.Pow(linear, 1f / 2.4f) - 0.055f;

        return (byte)Math.Clamp((int)MathF.Round(srgb * 255f), 0, 255);
    }
    
    /// <summary>
    /// Writes a linear RGB color to a byte buffer in BGRA32 format.
    /// </summary>
    /// <param name="buf">The byte array to write to (pixel buffer).</param>
    /// <param name="idx">The starting index in the buffer (must be 4-byte aligned).</param>
    /// <param name="linearRgb">The color in linear RGB space (components typically in [0, 1]).</param>
    /// <remarks>
    /// BGRA32 format stores pixels as: [Blue, Green, Red, Alpha] bytes.
    /// This is the native format for many graphics APIs and provides optimal performance.
    /// </remarks>
    private static void WriteBGRA(byte[] buf, int idx, Vector3 linearRgb)
    {
        buf[idx + 0] = ToSrgb8(linearRgb.Z); // Blue
        buf[idx + 1] = ToSrgb8(linearRgb.Y); // Green
        buf[idx + 2] = ToSrgb8(linearRgb.X); // Red
        buf[idx + 3] = 255;                   // Alpha
    }
    
    /// <summary>
    /// Renders a horizontal color gradient from red to green.
    /// Performs interpolation in linear RGB space and converts to sRGB for display.
    /// </summary>
    /// <remarks>
    /// This method demonstrates the correct way to interpolate colors:
    /// 1. Define color endpoints in linear RGB space
    /// 2. Perform linear interpolation in this perceptually uniform space
    /// 3. Convert the final linear RGB values to sRGB for display
    /// 
    /// Interpolating directly in sRGB space would produce incorrect, darker results
    /// because sRGB is not a linear color space.
    /// </remarks>
    private void DrawGradientLinearRgbThenToSrgb()
    {
        // Define gradient endpoints in Linear RGB space
        Vector3 leftColor  = new(1f, 0f, 0f); // Pure red
        Vector3 rightColor = new(0f, 1f, 0f); // Pure green

        // Get the stride for the BGRA32 format
        int stride =_writeableBitmap.BackBufferStride;
        
        // Allocate pixel buffer
        byte[] pixels = new byte[stride * Height];

        // Render the gradient
        for (int y = 0; y < Height; y++)
        {
            for (int x = 0; x < Width; x++)
            {
                // Calculate interpolation parameter (0.0 at left edge, 1.0 at right edge)
                float t = (float)x / (Width - 1);
                
                // Interpolate in linear RGB space
                Vector3 linearColor =Lerp(leftColor, rightColor,t);
                
                // Calculate buffer index for this pixel
                int pixelIndex = y * stride + x * 4;
                
                // Convert to sRGB and write to buffer
                WriteBGRA(pixels, pixelIndex, linearColor);
            }
            
        }

        // Write the pixel buffer to the bitmap
        _writeableBitmap.Lock();
        _writeableBitmap.WritePixels(new Int32Rect(0, 0, Width, Height), pixels, stride, 0);
        _writeableBitmap.AddDirtyRect(new Int32Rect(0, 0, Width, Height));
        _writeableBitmap.Unlock();
    }
}