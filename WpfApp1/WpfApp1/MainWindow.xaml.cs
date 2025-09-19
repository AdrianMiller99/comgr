using System.Numerics;
using System.Windows;
using System.Windows.Media;
using System.Windows.Media.Imaging;


namespace WpfApp1;

/// <summary>
/// Interaction logic for MainWindow.xaml
/// </summary>
public partial class MainWindow : Window
{
    const int W = 512, H = 256;
    private WriteableBitmap _wb;
    
    public MainWindow()
    {
        InitializeComponent();
        _wb = new WriteableBitmap(W, H, 96, 96, PixelFormats.Bgra32, null);
        Canvas.Source = _wb;
        DrawGradientLinearRgbThenToSrgb();
    }
    
    public static Vector3 Lerp(Vector3 a, Vector3 b, float t)
    {
        return a + (b - a) * t;
    }
    
    static byte ToSrgb8(float linear)
    {
        // 1) clip to [0..1]
        if (linear <= 0f) return 0;
        if (linear >= 1f) return 255;

        // 2) sRGB transfer function (gamma ≈ 2.4 piecewise is the standard)
        // Linear→sRGB
        float srgb = linear <= 0.0031308f
            ? 12.92f * linear
            : 1.055f * MathF.Pow(linear, 1f / 2.4f) - 0.055f;

        // 3) ×255
        int v = (int)MathF.Round(srgb * 255f);
        if (v < 0) v = 0; if (v > 255) v = 255;
        return (byte)v;
    }
    
    static void WriteBgra(byte[] buf, int idx, Vector3 linearRgb)
    {
        byte r = ToSrgb8(linearRgb.X);
        byte g = ToSrgb8(linearRgb.Y);
        byte b = ToSrgb8(linearRgb.Z);
        buf[idx + 0] = b;   // B
        buf[idx + 1] = g;   // G
        buf[idx + 2] = r;   // R
        buf[idx + 3] = 255; // A
    }
    
    void DrawGradientLinearRgbThenToSrgb()
    {
        // Endpoints in *Linear RGB* (not sRGB!)
        Vector3 left  = new(1f, 0f, 0f); // pure red
        Vector3 right = new(0f, 1f, 0f); // pure green

        int stride = _wb.BackBufferStride; // BGRA32 stride (see slides on stride/formats)
        byte[] pixels = new byte[stride * H];

        for (int y = 0; y < H; y++)
        {
            for (int x = 0; x < W; x++)
            {
                float t = (float)x / (W - 1);
                Vector3 linear = Lerp(left, right, t);  // interpolation in Linear RGB (lab requirement)
                int idx = y * stride + x * 4;
                WriteBgra(pixels, idx, linear);         // convert to sRGB 8-bit when writing
            }
        }

        _wb.Lock();
        _wb.WritePixels(new Int32Rect(0, 0, W, H), pixels, stride, 0);
        _wb.AddDirtyRect(new Int32Rect(0, 0, W, H));
        _wb.Unlock();
    }
}