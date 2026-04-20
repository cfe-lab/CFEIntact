import { Box, Text } from "@chakra-ui/react";

interface OutputPanelProps {
  title: string;
  content: string;
  colorScheme?: string;
}

const OutputPanel = ({ title, content, colorScheme = "gray" }: OutputPanelProps) => {
  const hasContent = content.trim().length > 0;
  
  return (
    <Box>
      <Text fontWeight="bold" mb={2}>{title}</Text>
      <Box
        p={4}
        bg={`${colorScheme}.50`}
        borderRadius="md"
        borderWidth={1}
        borderColor={`${colorScheme}.200`}
        maxH="300px"
        overflowY="auto"
      >
        {hasContent ? (
          <pre
            style={{
              margin: 0,
              whiteSpace: "pre-wrap",
              wordBreak: "break-word",
              fontFamily: "monospace",
              fontSize: "0.875rem",
            }}
          >
            {content}
          </pre>
        ) : (
          <Text color="gray.500" fontStyle="italic">
            (no output)
          </Text>
        )}
      </Box>
    </Box>
  );
};

export default OutputPanel;
